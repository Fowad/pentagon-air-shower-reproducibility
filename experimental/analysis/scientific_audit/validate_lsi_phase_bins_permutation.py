#!/usr/bin/env python3
"""Direct, checkpointed V12 permutation test of LS I orbital phase.

Ten disjoint phase-bin counts are generated in every randomized sky.  The two
published broad windows and the full-archive count are then exact sums of those
bins, retaining all correlations without an extra Monte Carlo search.  Exact-
circle block matrices are bit packed to keep the memory footprint small.  Every
trial uses the definitive V12 deterministic per-realization RNG convention.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd

from exact_permutation_moments import ROOT, V12_PATH, load_v12
from lsi_orbital_phase_exact import (
    ORBIT_DAYS,
    PHASE_EPOCH_MJD,
    SOURCE_DEC_DEG,
    SOURCE_RA_DEG,
)


STREAM = 61_399
LABELS = np.asarray(
    [f"phase_{index / 10:.1f}_{(index + 1) / 10:.1f}" for index in range(10)],
    dtype=str,
)

_WORKER_PREPARED = None
_WORKER_BASE_SEED = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def phase_bins(events: pd.DataFrame) -> np.ndarray:
    times = pd.to_datetime(events["utc_time"], utc=True, format="mixed")
    mjd = times.astype("int64").to_numpy(float) / 1.0e9 / 86400.0 + 40587.0
    phase = np.mod(mjd - PHASE_EPOCH_MJD, ORBIT_DAYS) / ORBIT_DAYS
    return np.minimum(np.floor(phase * 10).astype(np.int8), 9)


def packed_lookup(
    packed: np.ndarray, rows: np.ndarray, columns: np.ndarray
) -> np.ndarray:
    return ((packed[rows, columns >> 3] >> (columns & 7)) & 1).astype(bool)


def trial_chunk(task):
    """Evaluate a contiguous deterministic trial range in one worker."""
    start, stop = task
    if _WORKER_PREPARED is None or _WORKER_BASE_SEED is None:
        raise RuntimeError("Permutation worker was not initialized.")
    output = np.zeros((stop - start, len(LABELS)), dtype=np.int32)
    for output_row, trial in enumerate(range(start, stop)):
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [int(_WORKER_BASE_SEED), int(STREAM), int(trial)]
            )
        )
        counts = np.zeros(len(LABELS), dtype=np.int32)
        for packed, local_bins, rows in _WORKER_PREPARED:
            permutation = rng.permutation(len(local_bins))
            inside = packed_lookup(packed, rows, permutation)
            counts += np.bincount(
                local_bins[permutation][inside], minlength=len(LABELS)
            ).astype(np.int32)
        output[output_row] = counts
    return start, output


def prepare_blocks(v12, state, bins, radius_deg):
    prepared = []
    packed_bytes = 0
    for block_index, indices in enumerate(state.blocks):
        indices = np.asarray(indices, dtype=np.int64)
        local_bins = bins[indices]
        matrix = v12._exact_target_block_matrix(
            state,
            indices,
            SOURCE_RA_DEG,
            SOURCE_DEC_DEG,
            radius_deg,
        )
        packed = np.packbits(matrix, axis=1, bitorder="little")
        prepared.append(
            (packed, local_bins, np.arange(len(indices), dtype=np.int64))
        )
        packed_bytes += packed.nbytes
        if (block_index + 1) % 100 == 0:
            print(
                f"prepared through block {block_index + 1}/{len(state.blocks)}",
                flush=True,
            )
    print(
        f"prepared {len(prepared)} relevant blocks "
        f"({packed_bytes / 1024**2:.2f} MiB packed)",
        flush=True,
    )
    return prepared


def save_checkpoint(path, values, completed, trials, event_sha, state_sha, v12_sha):
    np.savez_compressed(
        path,
        values=values,
        completed=np.asarray(completed, dtype=np.int64),
        target_trials=np.asarray(trials, dtype=np.int64),
        labels=LABELS,
        event_sha256=np.asarray(event_sha),
        state_fingerprint=np.asarray(state_sha),
        v12_code_sha256=np.asarray(v12_sha),
        stream=np.asarray(STREAM, dtype=np.int64),
    )


def load_checkpoint(path, trials, event_sha, state_sha, v12_sha):
    if not path.exists():
        return np.zeros((trials, len(LABELS)), dtype=np.int32), 0
    with np.load(path, allow_pickle=False) as data:
        valid = (
            int(data["target_trials"]) == trials
            and str(data["event_sha256"]) == event_sha
            and str(data["state_fingerprint"]) == state_sha
            and str(data["v12_code_sha256"]) == v12_sha
            and int(data["stream"]) == STREAM
            and np.array_equal(data["labels"].astype(str), LABELS)
        )
        if not valid:
            raise RuntimeError(
                f"Checkpoint identity does not match this run: {path}"
            )
        values = np.asarray(data["values"], dtype=np.int32)
        completed = int(data["completed"])
    if values.shape != (trials, len(LABELS)) or not 0 <= completed <= trials:
        raise RuntimeError(f"Invalid checkpoint shape or completion count: {path}")
    print(f"resuming direct permutations at {completed}/{trials}", flush=True)
    return values, completed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--events",
        type=Path,
        default=ROOT / "audit" / "v12" / "reconstructed_events.csv.gz",
    )
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--checkpoint-every", type=int, default=5_000)
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, os.cpu_count() or 1),
        help="Parallel fork workers; use 1 for serial execution.",
    )
    parser.add_argument(
        "--trial-chunk",
        type=int,
        default=100,
        help="Contiguous trials returned by each worker task.",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=ROOT / "audit" / "results" / "lsi_phase_bins_100k_checkpoint.npz",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "audit" / "results" / "lsi_phase_bins_100k_summary.json",
    )
    args = parser.parse_args()
    if (
        args.trials <= 0
        or args.checkpoint_every <= 0
        or args.workers <= 0
        or args.trial_chunk <= 0
    ):
        raise ValueError("Trial, checkpoint, worker, and chunk counts must be positive.")

    v12 = load_v12()
    events = pd.read_csv(args.events)
    bins = phase_bins(events)
    config = v12.Config(
        data_dir=Path("."),
        output_dir=Path("."),
        event_cache=None,
        run_summary_cache=None,
        make_figures=False,
    )
    state = v12.randomization_state(
        events, config, np.zeros(len(events), dtype=np.int16)
    )
    event_sha = v12.event_fingerprint(events)
    state_sha = v12.randomization_state_fingerprint(state)
    v12_sha = sha256(V12_PATH)
    prepared = prepare_blocks(v12, state, bins, config.aperture_radius_deg)

    inside_observed = v12._exact_circle_membership(
        events["ra_deg"].to_numpy(float),
        events["dec_deg"].to_numpy(float),
        SOURCE_RA_DEG,
        SOURCE_DEC_DEG,
        config.aperture_radius_deg,
    )
    observed = np.bincount(
        bins[inside_observed], minlength=len(LABELS)
    ).astype(np.int32)

    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    values, completed = load_checkpoint(
        args.checkpoint, args.trials, event_sha, state_sha, v12_sha
    )
    tasks = [
        (start, min(start + args.trial_chunk, args.trials))
        for start in range(completed, args.trials, args.trial_chunk)
    ]
    global _WORKER_PREPARED, _WORKER_BASE_SEED
    _WORKER_PREPARED = prepared
    _WORKER_BASE_SEED = int(config.random_seed)

    last_checkpoint = completed
    if args.workers == 1:
        iterator = map(trial_chunk, tasks)
        pool = None
    else:
        if "fork" not in mp.get_all_start_methods():
            raise RuntimeError(
                "This optimized run requires the multiprocessing 'fork' start method. "
                "Use --workers 1 on platforms without it."
            )
        context = mp.get_context("fork")
        pool = context.Pool(processes=args.workers)
        iterator = pool.imap(trial_chunk, tasks, chunksize=1)

    try:
        for start, chunk_values in iterator:
            stop = start + len(chunk_values)
            values[start:stop] = chunk_values
            completed = stop
            if (
                completed - last_checkpoint >= args.checkpoint_every
                or completed == args.trials
            ):
                save_checkpoint(
                    args.checkpoint,
                    values,
                    completed,
                    args.trials,
                    event_sha,
                    state_sha,
                    v12_sha,
                )
                last_checkpoint = completed
                print(
                    f"direct permutations {completed}/{args.trials} "
                    f"with {args.workers} worker(s)",
                    flush=True,
                )
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    background = values.mean(axis=0)
    background_std = values.std(axis=0, ddof=1)
    observed_z = (observed - background) / background_std
    null_z = (values - background) / background_std
    observed_bin_max_z = float(np.max(observed_z))
    local_exceedances = np.count_nonzero(
        values >= observed[None, :], axis=0
    )
    local_p = (1 + local_exceedances) / (args.trials + 1)
    ten_bin_family_exceedances = int(
        np.count_nonzero(np.max(null_z, axis=1) >= observed_bin_max_z)
    )
    ten_bin_family_p = float(
        (1 + ten_bin_family_exceedances)
        / (args.trials + 1)
    )

    published_labels = np.asarray(
        [
            "published_25_100TeV_phase_0.3_0.6",
            "published_gt100TeV_phase_0.6_0.2",
        ],
        dtype=str,
    )
    published_observed = np.asarray(
        [observed[3:6].sum(), observed[[6, 7, 8, 9, 0, 1]].sum()],
        dtype=np.int32,
    )
    published_values = np.column_stack(
        [
            values[:, 3:6].sum(axis=1),
            values[:, [6, 7, 8, 9, 0, 1]].sum(axis=1),
        ]
    )
    published_background = published_values.mean(axis=0)
    published_std = published_values.std(axis=0, ddof=1)
    published_z = (
        published_observed - published_background
    ) / published_std
    published_null_z = (
        published_values - published_background
    ) / published_std
    published_local_exceedances = np.count_nonzero(
        published_values >= published_observed[None, :], axis=0
    )
    published_local_p = (1 + published_local_exceedances) / (args.trials + 1)
    published_max_z = float(np.max(published_z))
    published_family_exceedances = int(
        np.count_nonzero(
            np.max(published_null_z, axis=1) >= published_max_z
        )
    )
    published_family_p = float(
        (1 + published_family_exceedances) / (args.trials + 1)
    )

    full_observed = int(observed.sum())
    full_values = values.sum(axis=1)
    full_background = float(full_values.mean())
    full_std = float(full_values.std(ddof=1))
    full_z = float((full_observed - full_background) / full_std)
    full_exceedances = int(np.count_nonzero(full_values >= full_observed))
    full_p = float((1 + full_exceedances) / (args.trials + 1))
    result = {
        "release_id": "PENTAGON-PEV-FOLLOWUP-V12A-2026-08-20",
        "source": "LS I +61 303",
        "ra_deg": SOURCE_RA_DEG,
        "dec_deg": SOURCE_DEC_DEG,
        "phase_epoch_mjd": PHASE_EPOCH_MJD,
        "period_days": ORBIT_DAYS,
        "trials": args.trials,
        "workers": args.workers,
        "trial_chunk": args.trial_chunk,
        "randomization": "within-run four-hour time-label permutation",
        "aperture_radius_deg": config.aperture_radius_deg,
        "event_sha256": event_sha,
        "state_fingerprint": state_sha,
        "v12_code_sha256": v12_sha,
        "phase_bin_tests": [
            {
                "test": str(LABELS[index]),
                "N": int(observed[index]),
                "B": float(background[index]),
                "sB": float(background_std[index]),
                "Z": float(observed_z[index]),
                "exceedances": int(local_exceedances[index]),
                "one_sided_p": float(local_p[index]),
            }
            for index in range(len(LABELS))
        ],
        "ten_phase_bin_max_z": observed_bin_max_z,
        "ten_phase_bin_family_exceedances": ten_bin_family_exceedances,
        "ten_phase_bin_family_p": ten_bin_family_p,
        "published_window_tests": [
            {
                "test": str(published_labels[index]),
                "N": int(published_observed[index]),
                "B": float(published_background[index]),
                "sB": float(published_std[index]),
                "Z": float(published_z[index]),
                "exceedances": int(published_local_exceedances[index]),
                "one_sided_p": float(published_local_p[index]),
            }
            for index in range(len(published_labels))
        ],
        "two_published_window_max_z": published_max_z,
        "two_published_window_family_exceedances": published_family_exceedances,
        "two_published_window_family_p": published_family_p,
        "full_archive": {
            "N": full_observed,
            "B": full_background,
            "sB": full_std,
            "Z": full_z,
            "exceedances": full_exceedances,
            "one_sided_p": full_p,
        },
        "checkpoint": str(args.checkpoint),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
