#!/usr/bin/env python3
"""Exact first and second moments for source-centered V12 permutation counts.

For a fixed source and one V12 four-hour block, the randomized aperture count
is a sum over a uniformly random permutation of the block's time labels.  Its
mean and variance follow exactly from the corresponding binary membership
matrix.  Independent blocks add.  This avoids Monte Carlo noise when auditing
the full-archive source-centered Z values.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
V12_PATH = (
    ROOT
    / "audit"
    / "handoff"
    / "PENTAGON_V12_NEW_WORK_HANDOFF_PACKAGE"
    / "Pentagon_Array_Directional_Analysis_V12.py"
)


def load_v12():
    spec = importlib.util.spec_from_file_location("pentagon_v12", V12_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {V12_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def source_vectors(sources: pd.DataFrame) -> np.ndarray:
    ra = np.radians(sources["ra_deg"].to_numpy(float))
    dec = np.radians(sources["dec_deg"].to_numpy(float))
    return np.column_stack(
        (np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))
    )


def observed_counts(v12, events: pd.DataFrame, sources: pd.DataFrame) -> np.ndarray:
    return np.asarray(
        [
            np.count_nonzero(
                v12._exact_circle_membership(
                    events["ra_deg"].to_numpy(float),
                    events["dec_deg"].to_numpy(float),
                    row.ra_deg,
                    row.dec_deg,
                    8.0,
                )
            )
            for row in sources.itertuples(index=False)
        ],
        dtype=np.int64,
    )


def observed_counts_by_period(
    v12, events: pd.DataFrame, sources: pd.DataFrame, interval_id: np.ndarray
) -> np.ndarray:
    out = np.zeros((7, len(sources)), dtype=np.int64)
    for column, row in enumerate(sources.itertuples(index=False)):
        membership = v12._exact_circle_membership(
            events["ra_deg"].to_numpy(float),
            events["dec_deg"].to_numpy(float),
            row.ra_deg,
            row.dec_deg,
            8.0,
        )
        out[0, column] = np.count_nonzero(membership)
        for interval in range(6):
            out[interval + 1, column] = np.count_nonzero(
                membership & (interval_id == interval)
            )
    return out


def exact_moments(v12, state, sources: pd.DataFrame, chunk_size: int = 12):
    vectors = source_vectors(sources)
    source_count = len(sources)
    means = np.zeros(source_count, dtype=np.float64)
    variances = np.zeros(source_count, dtype=np.float64)
    cos_radius = math.cos(math.radians(8.0)) - 1.0e-15

    for block_index, indices in enumerate(state.blocks):
        indices = np.asarray(indices, dtype=np.int64)
        n = len(indices)
        hour_angle = state.hour_angle_deg[indices][:, None]
        dec_mean = np.broadcast_to(state.declination_mean_deg[indices][:, None], (n, n))
        lst = np.broadcast_to(state.observed_lst_deg[indices][None, :], (n, n))
        angles = [
            np.broadcast_to(angle[indices][None, :], (n, n))
            for angle in state.precession
        ]
        ra, dec = v12.mean_of_date_to_icrs_with_angles(
            (lst - hour_angle) % 360.0, dec_mean, *angles
        )
        ra = np.radians(ra.ravel())
        dec = np.radians(dec.ravel())
        sky_vectors = np.column_stack(
            (np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))
        )

        for start in range(0, source_count, chunk_size):
            stop = min(start + chunk_size, source_count)
            membership = (sky_vectors @ vectors[start:stop].T) >= cos_radius
            membership = membership.reshape(n, n, stop - start)
            total = membership.sum(axis=(0, 1), dtype=np.int64).astype(float)
            row_sum = membership.sum(axis=1, dtype=np.int64).astype(float)
            col_sum = membership.sum(axis=0, dtype=np.int64).astype(float)
            means[start:stop] += total / n
            if n > 1:
                centered_ss = (
                    total
                    - np.square(row_sum).sum(axis=0) / n
                    - np.square(col_sum).sum(axis=0) / n
                    + np.square(total) / (n * n)
                )
                variances[start:stop] += centered_ss / (n - 1)

        if block_index % 100 == 0 or block_index + 1 == len(state.blocks):
            print(f"blocks {block_index + 1}/{len(state.blocks)}", flush=True)

    return means, variances


def exact_period_moments(v12, state, sources: pd.DataFrame, chunk_size: int = 12):
    """Return exact moments for the full archive and six reassigned intervals."""
    vectors = source_vectors(sources)
    source_count = len(sources)
    means = np.zeros((7, source_count), dtype=np.float64)
    variances = np.zeros((7, source_count), dtype=np.float64)
    cos_radius = math.cos(math.radians(8.0)) - 1.0e-15

    for block_index, indices in enumerate(state.blocks):
        indices = np.asarray(indices, dtype=np.int64)
        n = len(indices)
        hour_angle = state.hour_angle_deg[indices][:, None]
        dec_mean = np.broadcast_to(state.declination_mean_deg[indices][:, None], (n, n))
        lst = np.broadcast_to(state.observed_lst_deg[indices][None, :], (n, n))
        angles = [
            np.broadcast_to(angle[indices][None, :], (n, n))
            for angle in state.precession
        ]
        ra, dec = v12.mean_of_date_to_icrs_with_angles(
            (lst - hour_angle) % 360.0, dec_mean, *angles
        )
        ra = np.radians(ra.ravel())
        dec = np.radians(dec.ravel())
        sky_vectors = np.column_stack(
            (np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))
        )
        period_masks = [np.ones(n, dtype=bool)] + [
            state.interval_id[indices] == interval for interval in range(6)
        ]

        for start in range(0, source_count, chunk_size):
            stop = min(start + chunk_size, source_count)
            membership = (sky_vectors @ vectors[start:stop].T) >= cos_radius
            membership = membership.reshape(n, n, stop - start)
            for period, column_mask in enumerate(period_masks):
                if not np.any(column_mask):
                    continue
                selected = membership & column_mask[None, :, None]
                total = selected.sum(axis=(0, 1), dtype=np.int64).astype(float)
                row_sum = selected.sum(axis=1, dtype=np.int64).astype(float)
                col_sum = selected.sum(axis=0, dtype=np.int64).astype(float)
                means[period, start:stop] += total / n
                if n > 1:
                    centered_ss = (
                        total
                        - np.square(row_sum).sum(axis=0) / n
                        - np.square(col_sum).sum(axis=0) / n
                        + np.square(total) / (n * n)
                    )
                    variances[period, start:stop] += centered_ss / (n - 1)

        if block_index % 100 == 0 or block_index + 1 == len(state.blocks):
            print(f"blocks {block_index + 1}/{len(state.blocks)}", flush=True)

    return means, variances


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument(
        "--events",
        type=Path,
        default=ROOT / "audit" / "v12" / "reconstructed_events.csv.gz",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=12)
    parser.add_argument("--include-intervals", action="store_true")
    args = parser.parse_args()

    v12 = load_v12()
    sources = pd.read_csv(args.sources)
    required = {"source", "ra_deg", "dec_deg"}
    if not required.issubset(sources.columns):
        raise ValueError(f"Source CSV needs columns {sorted(required)}")
    events = pd.read_csv(args.events)
    config = v12.Config(
        data_dir=Path("."),
        output_dir=Path("."),
        event_cache=None,
        run_summary_cache=None,
        make_figures=False,
    )
    if args.include_intervals:
        primary = json.loads(
            (ROOT / "audit" / "v12" / "results" / "primary_results_V12.json").read_text()
        )
        boundaries = pd.to_datetime(primary["boundaries"], utc=True, format="mixed")
        times = pd.to_datetime(events["utc_time"], utc=True, format="mixed")
        interval_id = np.searchsorted(
            boundaries.asi8[1:-1], times.astype("int64"), side="right"
        ).astype(np.int16)
        state = v12.randomization_state(events, config, interval_id)
        observed = observed_counts_by_period(v12, events, sources, interval_id)
        mean, variance = exact_period_moments(v12, state, sources, args.chunk_size)
        std = np.sqrt(np.maximum(variance, 0.0))
        frames = []
        for period in range(7):
            frame = sources.copy()
            frame.insert(0, "period", "full" if period == 0 else f"interval_{period}")
            frame["N"] = observed[period]
            frame["B_exact"] = mean[period]
            frame["sB_exact"] = std[period]
            frame["Z_exact_moments"] = np.divide(
                observed[period] - mean[period],
                std[period],
                out=np.full(len(frame), np.nan),
                where=std[period] > 0,
            )
            frames.append(frame)
        out = pd.concat(frames, ignore_index=True)
    else:
        state = v12.randomization_state(
            events, config, np.zeros(len(events), dtype=np.int16)
        )
        observed = observed_counts(v12, events, sources)
        mean, variance = exact_moments(v12, state, sources, args.chunk_size)
        std = np.sqrt(np.maximum(variance, 0.0))
        out = sources.copy()
        out["N"] = observed
        out["B_exact"] = mean
        out["sB_exact"] = std
        out["Z_exact_moments"] = np.divide(
            observed - mean, std, out=np.full(len(out), np.nan), where=std > 0
        )
    out["normal_one_sided_p"] = norm.sf(out["Z_exact_moments"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
