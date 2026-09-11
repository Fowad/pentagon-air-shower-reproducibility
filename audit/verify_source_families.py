#!/usr/bin/env python3
"""Independent recalculation of the manuscript's directed-source probabilities."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "recovered" / "scientific_audit_outputs" / "results"
OUT = ROOT / "fresh_evidence" / "source_family_probability_recheck.json"


def add_one(values: np.ndarray, observed: float) -> tuple[float, int]:
    exceed = int(np.count_nonzero(np.asarray(values) >= observed))
    return (exceed + 1.0) / (len(values) + 1.0), exceed


def standardized(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = values.mean(axis=0)
    std = values.std(axis=0, ddof=1)
    z = (values - mean) / std
    return z, mean, std


def gaussian_max(
    correlation: np.ndarray,
    threshold: float,
    *,
    seed: int,
    trials: int = 1_000_000,
) -> tuple[float, int, float]:
    correlation = (correlation + correlation.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    factor = eigenvectors @ np.diag(np.sqrt(np.clip(eigenvalues, 0.0, None)))
    rng = np.random.default_rng(seed)
    exceed = 0
    for start in range(0, trials, 20_000):
        count = min(20_000, trials - start)
        draws = rng.standard_normal((count, correlation.shape[0])) @ factor.T
        exceed += int(np.count_nonzero(draws.max(axis=1) >= threshold))
    return (exceed + 1.0) / (trials + 1.0), exceed, float(eigenvalues.min())


def main() -> None:
    expected_event = "f7f370ad96d69ef9a793472490b4f10f4deff4e6c62cd2211527f1c65e245a6b"
    expected_code = "853c690442f3f956661bb41707192cbb1a9438a6f7be17e13f85c48bb1ea3435"
    report: dict = {"checks": {}, "direct_empirical": {}, "gaussian_covariance": {}}

    # LS I interval-1 local probability.
    with np.load(RESULTS / "lsi_interval1_100k_checkpoint.npz", allow_pickle=False) as d:
        values = d["values"].copy()
        identity = json.loads(str(d["identity_json"]))
        report["checks"]["lsi_interval_checkpoint_complete"] = int(d["completed"]) == 100_000
        report["checks"]["lsi_interval_identity"] = (
            identity["event_sha256"] == expected_event and identity["v12_code_sha256"] == expected_code
        )
    p, exceed = add_one(values, 854)
    report["direct_empirical"]["lsi_interval1_local"] = {
        "observed_count": 854,
        "trials": len(values),
        "exceedances": exceed,
        "p": p,
    }

    # Full archive, ten phase bins, and the two predeclared broad phase windows.
    phase_summary = json.loads((RESULTS / "lsi_phase_bins_100k_summary.json").read_text())
    with np.load(RESULTS / "lsi_phase_bins_100k_checkpoint.npz", allow_pickle=False) as d:
        phase_values = d["values"].copy()
        report["checks"]["phase_checkpoint_complete"] = (
            int(d["completed"]) == int(d["target_trials"]) == 100_000
        )
        report["checks"]["phase_checkpoint_identity"] = (
            str(d["event_sha256"]) == expected_event and str(d["v12_code_sha256"]) == expected_code
        )
    phase_z, phase_mean, phase_std = standardized(phase_values)
    observed_phase = np.asarray([row["N"] for row in phase_summary["phase_bin_tests"]], float)
    observed_phase_z = (observed_phase - phase_mean) / phase_std
    ten_p, ten_exceed = add_one(phase_z.max(axis=1), observed_phase_z.max())
    full_null = phase_values.sum(axis=1)
    full_p, full_exceed = add_one(full_null, phase_summary["full_archive"]["N"])
    broad_null = np.column_stack(
        (
            phase_values[:, 3:6].sum(axis=1),
            phase_values[:, [6, 7, 8, 9, 0, 1]].sum(axis=1),
        )
    )
    broad_z, broad_mean, broad_std = standardized(broad_null)
    broad_observed = np.asarray([row["N"] for row in phase_summary["published_window_tests"]], float)
    broad_observed_z = (broad_observed - broad_mean) / broad_std
    broad_p, broad_exceed = add_one(broad_z.max(axis=1), broad_observed_z.max())
    report["direct_empirical"].update(
        {
            "lsi_full_archive_local": {"p": full_p, "exceedances": full_exceed},
            "lsi_ten_phase_family": {"p": ten_p, "exceedances": ten_exceed},
            "lsi_two_published_phase_windows": {"p": broad_p, "exceedances": broad_exceed},
        }
    )

    # Six intervals x ten phase cells and the 12 interval/published-window tests.
    grid_summary = json.loads((RESULTS / "lsi_interval_phase_grid_100k_summary.json").read_text())
    with np.load(RESULTS / "lsi_interval_phase_grid_100k_checkpoint.npz", allow_pickle=False) as d:
        grid_values = d["values"].copy()
        report["checks"]["interval_phase_checkpoint_complete"] = (
            int(d["completed"]) == int(d["target_trials"]) == 100_000
        )
        report["checks"]["interval_phase_checkpoint_identity"] = (
            str(d["event_sha256"]) == expected_event and str(d["v12_code_sha256"]) == expected_code
        )
    grid_z, grid_mean, grid_std = standardized(grid_values)
    observed_grid = np.asarray(
        [row["N"] for row in grid_summary["sixty_interval_phase_cell_tests"]], float
    )
    observed_grid_z = (observed_grid - grid_mean) / grid_std
    grid_p, grid_exceed = add_one(grid_z.max(axis=1), observed_grid_z.max())
    grid_cube = grid_values.reshape(len(grid_values), 6, 10)
    published_interval_null = np.stack(
        (
            grid_cube[:, :, 3:6].sum(axis=2),
            grid_cube[:, :, [6, 7, 8, 9, 0, 1]].sum(axis=2),
        ),
        axis=2,
    ).reshape(len(grid_values), 12)
    pub_z, pub_mean, pub_std = standardized(published_interval_null)
    observed_pub = np.asarray(
        [row["N"] for row in grid_summary["interval_published_window_tests"]], float
    )
    observed_pub_z = (observed_pub - pub_mean) / pub_std
    pub_p, pub_exceed = add_one(pub_z.max(axis=1), observed_pub_z.max())
    interval_null = grid_cube.sum(axis=2)
    interval_z, interval_mean, interval_std = standardized(interval_null)
    observed_interval = observed_grid.reshape(6, 10).sum(axis=1)
    observed_interval_z = (observed_interval - interval_mean) / interval_std
    source_six_p, source_six_exceed = add_one(
        interval_z.max(axis=1), observed_interval_z.max()
    )
    report["direct_empirical"].update(
        {
            "lsi_sixty_interval_phase_cells": {"p": grid_p, "exceedances": grid_exceed},
            "lsi_twelve_interval_published_windows": {"p": pub_p, "exceedances": pub_exceed},
            "lsi_source_specific_six_intervals": {
                "p": source_six_p,
                "exceedances": source_six_exceed,
                "observed_max_z": float(observed_interval_z.max()),
            },
        }
    )

    # Recalibrate the two correlated Gaussian source families with new RNG streams.
    with np.load(RESULTS / "lhaaso_37_exact_full_covariance.npz", allow_pickle=False) as d:
        obs = d["observed"].astype(float)
        mean = d["mean"].copy()
        cov = d["covariance"].copy()
    std = np.sqrt(np.diag(cov))
    z = (obs - mean) / std
    corr = cov / np.outer(std, std)
    lhaaso_p, lhaaso_exceed, lhaaso_eig = gaussian_max(
        corr, float(z.max()), seed=2026090201
    )
    equal_z = float(np.ones(len(z)) @ z / np.sqrt(np.ones(len(z)) @ corr @ np.ones(len(z))))
    report["gaussian_covariance"]["lhaaso_37_full"] = {
        "observed_max_z": float(z.max()),
        "independent_trials": 1_000_000,
        "independent_p": lhaaso_p,
        "independent_exceedances": lhaaso_exceed,
        "minimum_correlation_eigenvalue": lhaaso_eig,
        "equal_source_stack_z": equal_z,
        "equal_source_stack_one_sided_p": float(norm.sf(equal_z)),
    }
    source_table = pd.read_csv(RESULTS.parent / "lhaaso_37_sources.csv")
    with np.load(RESULTS / "lhaaso_37_exact_full_covariance.npz", allow_pickle=False) as order_data:
        source_order = [str(item) for item in order_data["source"]]
    if source_table.source.astype(str).tolist() != source_order:
        raise RuntimeError("LHAASO source order differs between covariance and catalog")
    flux100 = source_table.flux_norm_50tev.to_numpy(float) * np.power(
        2.0, -source_table.spectral_index.to_numpy(float)
    )
    flux100_z = float(flux100 @ z / np.sqrt(flux100 @ corr @ flux100))
    report["gaussian_covariance"]["lhaaso_37_full"].update(
        {
            "flux100_weighted_stack_z": flux100_z,
            "flux100_weighted_stack_one_sided_p": float(norm.sf(flux100_z)),
        }
    )

    with np.load(RESULTS / "p300_exact_family_covariance.npz", allow_pickle=False) as d:
        source_count = len(d["source"])
        obs = d["observed"].astype(float)
        mean = d["mean"].copy()
        cov = d["covariance"].copy()
    std = np.sqrt(np.diag(cov))
    z = (obs - mean) / std
    corr = cov / np.outer(std, std)
    interval_p, interval_exceed, interval_eig = gaussian_max(
        corr, float(z.max()), seed=2026090202
    )
    aggregate = np.zeros((source_count, 6 * source_count))
    for interval in range(6):
        aggregate[:, interval * source_count : (interval + 1) * source_count] = np.eye(source_count)
    full_obs = aggregate @ obs
    full_mean = aggregate @ mean
    full_cov = aggregate @ cov @ aggregate.T
    full_std = np.sqrt(np.diag(full_cov))
    full_z = (full_obs - full_mean) / full_std
    full_corr = full_cov / np.outer(full_std, full_std)
    full_p, full_exceed, full_eig = gaussian_max(
        full_corr, float(full_z.max()), seed=2026090203
    )
    report["gaussian_covariance"]["uhe_17_by_6"] = {
        "observed_max_z": float(z.max()),
        "independent_trials": 1_000_000,
        "independent_p": interval_p,
        "independent_exceedances": interval_exceed,
        "minimum_correlation_eigenvalue": interval_eig,
    }
    report["gaussian_covariance"]["uhe_17_full"] = {
        "observed_max_z": float(full_z.max()),
        "independent_trials": 1_000_000,
        "independent_p": full_p,
        "independent_exceedances": full_exceed,
        "minimum_correlation_eigenvalue": full_eig,
    }

    # Two predeclared 17-source stacks across the six equal-live intervals.
    equal_transform = np.zeros((6, 6 * source_count), dtype=float)
    count_transform = np.zeros((6, 6 * source_count), dtype=float)
    for interval in range(6):
        block = slice(interval * source_count, (interval + 1) * source_count)
        weights = np.ones(source_count)
        denominator = np.sqrt(weights @ corr[block, block] @ weights)
        equal_transform[interval, block] = weights / denominator
        count_transform[interval, block] = weights
    equal_stack_z = equal_transform @ z
    equal_stack_corr = equal_transform @ corr @ equal_transform.T
    equal_stack_p, equal_stack_exceed, equal_stack_eig = gaussian_max(
        equal_stack_corr, float(equal_stack_z.max()), seed=2026090204
    )

    count_observed = count_transform @ obs
    count_mean = count_transform @ mean
    count_cov = count_transform @ cov @ count_transform.T
    count_std = np.sqrt(np.diag(count_cov))
    count_stack_z = (count_observed - count_mean) / count_std
    count_stack_corr = count_cov / np.outer(count_std, count_std)
    count_stack_p, count_stack_exceed, count_stack_eig = gaussian_max(
        count_stack_corr, float(count_stack_z.max()), seed=2026090205
    )
    report["gaussian_covariance"]["uhe_17_interval_stacks"] = {
        "equal_source_stack_z": equal_stack_z.tolist(),
        "equal_source_observed_max_z": float(equal_stack_z.max()),
        "equal_source_independent_trials": 1_000_000,
        "equal_source_independent_exceedances": equal_stack_exceed,
        "equal_source_independent_six_interval_p": equal_stack_p,
        "equal_source_minimum_correlation_eigenvalue": equal_stack_eig,
        "aggregate_count_stack_z": count_stack_z.tolist(),
        "aggregate_count_observed_max_z": float(count_stack_z.max()),
        "aggregate_count_independent_trials": 1_000_000,
        "aggregate_count_independent_exceedances": count_stack_exceed,
        "aggregate_count_independent_six_interval_p": count_stack_p,
        "aggregate_count_minimum_correlation_eigenvalue": count_stack_eig,
    }

    expected = {
        "lsi_interval1_local": 0.0015799842001579985,
        "lsi_full_archive_local": 0.02974970250297497,
        "lsi_ten_phase_family": 0.0930590694093059,
        "lsi_two_published_phase_windows": 0.031949680503194966,
        "lsi_sixty_interval_phase_cells": 0.11538884611153888,
        "lsi_twelve_interval_published_windows": 0.07383926160738392,
        "lsi_source_specific_six_intervals": 0.00974990250097499,
    }
    for key, value in expected.items():
        report["checks"][f"{key}_exact"] = abs(report["direct_empirical"][key]["p"] - value) < 1e-15
    report["checks"]["lhaaso_gaussian_agrees"] = abs(lhaaso_p - 0.7178658564268288) < 0.003
    report["checks"]["uhe_interval_gaussian_agrees"] = abs(interval_p - 0.13833137233372553) < 0.003
    report["checks"]["uhe_full_gaussian_agrees"] = abs(full_p - 0.3557955288408942) < 0.003
    report["checks"]["lhaaso_flux100_stack_exact"] = abs(
        report["gaussian_covariance"]["lhaaso_37_full"]["flux100_weighted_stack_one_sided_p"]
        - 0.7054416718753715
    ) < 1e-12
    report["checks"]["uhe_equal_source_stack_gaussian_agrees"] = abs(
        equal_stack_p - 0.13076117384776523
    ) < 0.003
    report["checks"]["uhe_aggregate_count_stack_gaussian_agrees"] = abs(
        count_stack_p - 0.15131976973604605
    ) < 0.003
    report["all_checks_pass"] = all(report["checks"].values())
    if not report["all_checks_pass"]:
        raise RuntimeError([key for key, value in report["checks"].items() if not value])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
