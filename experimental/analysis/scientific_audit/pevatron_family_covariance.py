#!/usr/bin/env python3
"""Correlated exact-circle PeVatron-family calibration from V12 permutations.

The covariance of two statistics formed from a common uniform random
permutation has a closed form.  Summing those blockwise covariances gives the
joint Gaussian limit for the externally fixed source-by-interval family.  The
limit is exceptionally well motivated here because each count sums over 1,279
independent four-hour blocks.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from exact_permutation_moments import ROOT, load_v12, source_vectors


def angular_separation_matrix(sources: pd.DataFrame) -> np.ndarray:
    vectors = source_vectors(sources)
    return np.degrees(
        np.arccos(np.clip(vectors @ vectors.T, -1.0, 1.0))
    )


def observed_by_interval(v12, events, sources, interval_id):
    out = np.zeros((6, len(sources)), dtype=np.int64)
    for source_index, row in enumerate(sources.itertuples(index=False)):
        inside = v12._exact_circle_membership(
            events["ra_deg"].to_numpy(float),
            events["dec_deg"].to_numpy(float),
            row.ra_deg,
            row.dec_deg,
            8.0,
        )
        for interval in range(6):
            out[interval, source_index] = np.count_nonzero(
                inside & (interval_id == interval)
            )
    return out


def blockwise_joint_moments(v12, state, sources):
    m = len(sources)
    q = 6 * m
    vectors = source_vectors(sources)
    separations = angular_separation_matrix(sources)
    overlapping_pairs = [
        (a, b)
        for a in range(m)
        for b in range(a, m)
        if a == b or separations[a, b] < 16.0 + 1.0e-12
    ]
    mean = np.zeros(q, dtype=float)
    covariance = np.zeros((q, q), dtype=float)
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
        membership = ((sky_vectors @ vectors.T) >= cos_radius).reshape(n, n, m)

        statistics = np.zeros((n, n, 6, m), dtype=bool)
        for interval in np.unique(state.interval_id[indices]):
            mask = state.interval_id[indices] == interval
            statistics[:, :, int(interval), :] = membership & mask[None, :, None]
        statistics = statistics.reshape(n, n, q)
        total = statistics.sum(axis=(0, 1), dtype=np.int64).astype(float)
        row_sum = statistics.sum(axis=1, dtype=np.int64).astype(float)
        col_sum = statistics.sum(axis=0, dtype=np.int64).astype(float)

        overlap = np.zeros((q, q), dtype=float)
        present_intervals = np.unique(state.interval_id[indices])
        for interval in present_intervals:
            base = int(interval) * m
            for a, b in overlapping_pairs:
                value = float(
                    np.count_nonzero(
                        statistics[:, :, base + a] & statistics[:, :, base + b]
                    )
                )
                overlap[base + a, base + b] = value
                overlap[base + b, base + a] = value

        mean += total / n
        if n > 1:
            residual_cross = (
                overlap
                - (row_sum.T @ row_sum) / n
                - (col_sum.T @ col_sum) / n
                + np.outer(total, total) / (n * n)
            )
            covariance += residual_cross / (n - 1)

        if block_index % 100 == 0 or block_index + 1 == len(state.blocks):
            print(f"blocks {block_index + 1}/{len(state.blocks)}", flush=True)

    return mean, covariance


def gaussian_max_probability(correlation, threshold, trials=5_000_000, seed=271828):
    eigenvalues, eigenvectors = np.linalg.eigh((correlation + correlation.T) / 2)
    clipped = np.maximum(eigenvalues, 0.0)
    factor = eigenvectors @ np.diag(np.sqrt(clipped))
    rng = np.random.default_rng(seed)
    exceed = 0
    completed = 0
    chunk = 50_000
    while completed < trials:
        count = min(chunk, trials - completed)
        values = rng.standard_normal((count, correlation.shape[0])) @ factor.T
        exceed += int(np.count_nonzero(values.max(axis=1) >= threshold))
        completed += count
    return (1 + exceed) / (trials + 1), exceed, eigenvalues.min()


def linear_stack(observed, mean, covariance, weights):
    weights = np.asarray(weights, dtype=float)
    numerator = float(weights @ (observed - mean))
    standard_deviation = float(np.sqrt(weights @ covariance @ weights))
    return numerator / standard_deviation


def six_interval_stack_family(observed, mean, covariance, weights, seed):
    """Calibrate the maximum of six correlated linear interval stacks."""
    weights = np.asarray(weights, dtype=float)
    stacked_observed = weights @ observed
    stacked_mean = weights @ mean
    stacked_covariance = weights @ covariance @ weights.T
    stacked_std = np.sqrt(np.diag(stacked_covariance))
    stacked_z = (stacked_observed - stacked_mean) / stacked_std
    stacked_correlation = stacked_covariance / np.outer(stacked_std, stacked_std)
    threshold = float(np.max(stacked_z))
    p_value, exceedances, minimum_eigenvalue = gaussian_max_probability(
        stacked_correlation, threshold, seed=seed
    )
    return stacked_z, threshold, p_value, exceedances, minimum_eigenvalue


def main():
    v12 = load_v12()
    source_table = pd.read_csv(ROOT / "audit" / "pevatron_sources_working.csv")
    sources = source_table[source_table["p300_core"].eq("yes")].reset_index(drop=True)
    events = pd.read_csv(ROOT / "audit" / "v12" / "reconstructed_events.csv.gz")
    primary = json.loads(
        (ROOT / "audit" / "v12" / "results" / "primary_results_V12.json").read_text()
    )
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
    observed = observed_by_interval(v12, events, sources, interval_id).reshape(-1)
    mean, covariance = blockwise_joint_moments(v12, state, sources)
    std = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    z = np.divide(observed - mean, std, out=np.full_like(mean, np.nan), where=std > 0)
    correlation = covariance / np.outer(std, std)
    correlation = np.nan_to_num(correlation)
    np.fill_diagonal(correlation, 1.0)

    interval_threshold = float(np.nanmax(z))
    interval_p, interval_exceed, min_eigenvalue = gaussian_max_probability(
        correlation, interval_threshold
    )

    aggregate = np.zeros((len(sources), 6 * len(sources)), dtype=float)
    for interval in range(6):
        aggregate[:, interval * len(sources) : (interval + 1) * len(sources)] = np.eye(
            len(sources)
        )
    observed_full = aggregate @ observed
    mean_full = aggregate @ mean
    covariance_full = aggregate @ covariance @ aggregate.T
    std_full = np.sqrt(np.diag(covariance_full))
    z_full = (observed_full - mean_full) / std_full
    correlation_full = covariance_full / np.outer(std_full, std_full)
    full_threshold = float(np.max(z_full))
    full_p, full_exceed, full_min_eigenvalue = gaussian_max_probability(
        correlation_full, full_threshold, seed=314159
    )

    full_aggregate_count_z = linear_stack(
        observed, mean, covariance, np.ones(6 * len(sources))
    )
    full_equal_source_weights = aggregate.T @ (1.0 / std_full)
    full_equal_source_z = linear_stack(
        observed, mean, covariance, full_equal_source_weights
    )

    interval_aggregate_weights = np.zeros((6, 6 * len(sources)), dtype=float)
    interval_equal_weights = np.zeros_like(interval_aggregate_weights)
    diagonal_std = np.sqrt(np.diag(covariance))
    for interval in range(6):
        sl = slice(interval * len(sources), (interval + 1) * len(sources))
        interval_aggregate_weights[interval, sl] = 1.0
        interval_equal_weights[interval, sl] = 1.0 / diagonal_std[sl]
    (
        interval_aggregate_z,
        interval_aggregate_threshold,
        interval_aggregate_p,
        interval_aggregate_exceedances,
        interval_aggregate_min_eigenvalue,
    ) = six_interval_stack_family(
        observed,
        mean,
        covariance,
        interval_aggregate_weights,
        seed=161803,
    )
    (
        interval_equal_z,
        interval_equal_threshold,
        interval_equal_p,
        interval_equal_exceedances,
        interval_equal_min_eigenvalue,
    ) = six_interval_stack_family(
        observed,
        mean,
        covariance,
        interval_equal_weights,
        seed=141421,
    )

    rows = []
    for interval in range(6):
        for source_index, source in sources.iterrows():
            index = interval * len(sources) + source_index
            rows.append(
                {
                    "period": f"interval_{interval + 1}",
                    "source": source["source"],
                    "N": int(observed[index]),
                    "B_exact": mean[index],
                    "sB_exact": std[index],
                    "Z_exact": z[index],
                }
            )
    for source_index, source in sources.iterrows():
        rows.append(
            {
                "period": "full",
                "source": source["source"],
                "N": int(observed_full[source_index]),
                "B_exact": mean_full[source_index],
                "sB_exact": std_full[source_index],
                "Z_exact": z_full[source_index],
            }
        )
    table = pd.DataFrame(rows)
    result_dir = ROOT / "audit" / "results"
    result_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(result_dir / "p300_exact_family_statistics.csv", index=False)
    np.savez_compressed(
        result_dir / "p300_exact_family_covariance.npz",
        source=np.asarray(sources["source"], dtype=str),
        mean=mean,
        covariance=covariance,
        observed=observed,
        z=z,
        interval_threshold=interval_threshold,
        interval_family_p=interval_p,
        full_threshold=full_threshold,
        full_family_p=full_p,
    )
    summary = {
        "source_count": len(sources),
        "interval_test_count": 6 * len(sources),
        "interval_observed_max_z": interval_threshold,
        "interval_gaussian_family_p": interval_p,
        "interval_gaussian_exceedances": interval_exceed,
        "interval_correlation_min_eigenvalue_before_clipping": min_eigenvalue,
        "full_observed_max_z": full_threshold,
        "full_gaussian_family_p": full_p,
        "full_gaussian_exceedances": full_exceed,
        "full_correlation_min_eigenvalue_before_clipping": full_min_eigenvalue,
        "full_aggregate_count_stack_z": full_aggregate_count_z,
        "full_equal_source_stack_z": full_equal_source_z,
        "interval_aggregate_count_stack_z": interval_aggregate_z.tolist(),
        "interval_aggregate_count_max_z": interval_aggregate_threshold,
        "interval_aggregate_count_six_interval_family_p": interval_aggregate_p,
        "interval_aggregate_count_exceedances": interval_aggregate_exceedances,
        "interval_aggregate_count_min_eigenvalue": interval_aggregate_min_eigenvalue,
        "interval_equal_source_stack_z": interval_equal_z.tolist(),
        "interval_equal_source_max_z": interval_equal_threshold,
        "interval_equal_source_six_interval_family_p": interval_equal_p,
        "interval_equal_source_exceedances": interval_equal_exceedances,
        "interval_equal_source_min_eigenvalue": interval_equal_min_eigenvalue,
        "gaussian_trials": 5_000_000,
    }
    (result_dir / "p300_exact_family_summary_v12a.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(json.dumps(summary, indent=2))
    print(table.sort_values("Z_exact", ascending=False).head(15).to_string(index=False))


if __name__ == "__main__":
    main()
