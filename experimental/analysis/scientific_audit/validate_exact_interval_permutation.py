#!/usr/bin/env python3
"""Checkpointed direct V12 permutation test for one exact source/interval."""

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


_PREPARED = None
_BASE_SEED = None
_STREAM = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def packed_lookup(packed, rows, columns):
    return ((packed[rows, columns >> 3] >> (columns & 7)) & 1).astype(bool)


def trial_chunk(task):
    start, stop = task
    if _PREPARED is None:
        raise RuntimeError("Worker state was not initialized.")
    output = np.zeros(stop - start, dtype=np.int32)
    for destination, trial in enumerate(range(start, stop)):
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [int(_BASE_SEED), int(_STREAM), int(trial)]
            )
        )
        count = 0
        for packed, target_labels, rows in _PREPARED:
            permutation = rng.permutation(len(target_labels))
            inside = packed_lookup(packed, rows, permutation)
            count += int(np.count_nonzero(inside & target_labels[permutation]))
        output[destination] = count
    return start, output


def save_checkpoint(path, values, completed, identity):
    np.savez_compressed(
        path,
        values=values,
        completed=np.asarray(completed, dtype=np.int64),
        identity_json=np.asarray(json.dumps(identity, sort_keys=True)),
    )


def load_checkpoint(path, trials, identity):
    if not path.exists():
        return np.zeros(trials, dtype=np.int32), 0
    with np.load(path, allow_pickle=False) as data:
        if str(data["identity_json"]) != json.dumps(identity, sort_keys=True):
            raise RuntimeError(f"Checkpoint identity mismatch: {path}")
        values = np.asarray(data["values"], dtype=np.int32)
        completed = int(data["completed"])
    if values.shape != (trials,) or not 0 <= completed <= trials:
        raise RuntimeError(f"Invalid checkpoint dimensions: {path}")
    print(f"resuming at {completed}/{trials}", flush=True)
    return values, completed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--ra-deg", type=float, required=True)
    parser.add_argument("--dec-deg", type=float, required=True)
    parser.add_argument("--interval", type=int, choices=range(1, 7), required=True)
    parser.add_argument("--stream", type=int, required=True)
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--checkpoint-every", type=int, default=5_000)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--trial-chunk", type=int, default=200)
    parser.add_argument(
        "--events",
        type=Path,
        default=ROOT / "audit" / "v12" / "reconstructed_events.csv.gz",
    )
    parser.add_argument(
        "--primary-results",
        type=Path,
        default=ROOT / "audit" / "v12" / "results" / "primary_results_V12.json",
    )
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    v12 = load_v12()
    events = pd.read_csv(args.events)
    primary = json.loads(args.primary_results.read_text())
    boundaries = pd.to_datetime(primary["boundaries"], utc=True, format="mixed")
    times = pd.to_datetime(events["utc_time"], utc=True, format="mixed")
    interval_id = np.searchsorted(
        boundaries.asi8[1:-1], times.astype("int64"), side="right"
    ).astype(np.int16)
    config = v12.Config(
        data_dir=Path("."),
        output_dir=Path("."),
        event_cache=None,
        run_summary_cache=None,
        make_figures=False,
    )
    state = v12.randomization_state(events, config, interval_id)
    target_interval = args.interval - 1
    prepared = []
    packed_bytes = 0
    for block_index, indices in enumerate(state.blocks):
        indices = np.asarray(indices, dtype=np.int64)
        target_labels = state.interval_id[indices] == target_interval
        if not np.any(target_labels):
            continue
        matrix = v12._exact_target_block_matrix(
            state,
            indices,
            args.ra_deg,
            args.dec_deg,
            config.aperture_radius_deg,
        )
        packed = np.packbits(matrix, axis=1, bitorder="little")
        prepared.append(
            (packed, target_labels, np.arange(len(indices), dtype=np.int64))
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

    observed_mask = interval_id == target_interval
    observed = int(
        np.count_nonzero(
            v12._exact_circle_membership(
                events.loc[observed_mask, "ra_deg"].to_numpy(float),
                events.loc[observed_mask, "dec_deg"].to_numpy(float),
                args.ra_deg,
                args.dec_deg,
                config.aperture_radius_deg,
            )
        )
    )
    identity = {
        "release_id": "PENTAGON-PEV-FOLLOWUP-V12A-2026-08-20",
        "source": args.source,
        "ra_deg": args.ra_deg,
        "dec_deg": args.dec_deg,
        "interval": args.interval,
        "stream": args.stream,
        "trials": args.trials,
        "event_sha256": v12.event_fingerprint(events),
        "state_fingerprint": v12.randomization_state_fingerprint(state),
        "v12_code_sha256": sha256(V12_PATH),
    }
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    values, completed = load_checkpoint(args.checkpoint, args.trials, identity)

    global _PREPARED, _BASE_SEED, _STREAM
    _PREPARED = prepared
    _BASE_SEED = int(config.random_seed)
    _STREAM = int(args.stream)
    tasks = [
        (start, min(start + args.trial_chunk, args.trials))
        for start in range(completed, args.trials, args.trial_chunk)
    ]
    last_checkpoint = completed
    if args.workers == 1:
        iterator = map(trial_chunk, tasks)
        pool = None
    else:
        context = mp.get_context("fork")
        pool = context.Pool(args.workers)
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
                save_checkpoint(args.checkpoint, values, completed, identity)
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

    background = float(values.mean())
    background_std = float(values.std(ddof=1))
    z_value = float((observed - background) / background_std)
    exceedances = int(np.count_nonzero(values >= observed))
    p_value = float((1 + exceedances) / (args.trials + 1))
    result = {
        **identity,
        "aperture_radius_deg": config.aperture_radius_deg,
        "randomization": "within-run four-hour time/interval-label permutation",
        "workers": args.workers,
        "N": observed,
        "B": background,
        "sB": background_std,
        "Z": z_value,
        "exceedances": exceedances,
        "one_sided_p": p_value,
        "checkpoint": str(args.checkpoint),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
