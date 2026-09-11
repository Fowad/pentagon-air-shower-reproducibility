#!/usr/bin/env python3
"""Locked full-archive UHE Galactic-plane template family.

Implements Test 2 in UHE_ANALYSIS_PREREGISTRATION.md.  Run only after the
published LHAASO longitude boundaries have been independently verified.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
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
EVENTS_PATH = ROOT / "audit" / "v12" / "reconstructed_events.csv.gz"
OUTPUT_DIR = ROOT / "investigation" / "results"
EXPECTED_EVENTS_SHA256 = "2ca83067b8b14ee0c1a82eb8b4c51eef63c4b1d2b6e0526c9ddaa7ef6f8ebb40"
TEMPLATE_LABELS = (
    "all_longitudes_abs_b_lt_5deg",
    "all_longitudes_abs_b_lt_10deg",
    "lhaaso_inner_15_lt_l_lt_125_abs_b_lt_5deg",
    "lhaaso_outer_125_lt_l_lt_235_abs_b_lt_5deg",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_v12():
    spec = importlib.util.spec_from_file_location("pentagon_v12_uhe_plane", V12_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {V12_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def template_membership(longitude_deg: np.ndarray, latitude_deg: np.ndarray) -> np.ndarray:
    longitude = np.asarray(longitude_deg, dtype=float) % 360.0
    absolute_latitude = np.abs(np.asarray(latitude_deg, dtype=float))
    return np.column_stack(
        (
            absolute_latitude < 5.0,
            absolute_latitude < 10.0,
            (longitude >= 15.0) & (longitude < 125.0) & (absolute_latitude < 5.0),
            (longitude >= 125.0) & (longitude < 235.0) & (absolute_latitude < 5.0),
        )
    )


def exact_joint_moments(v12, state) -> tuple[np.ndarray, np.ndarray]:
    template_count = len(TEMPLATE_LABELS)
    mean = np.zeros(template_count, dtype=float)
    covariance = np.zeros((template_count, template_count), dtype=float)

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
        longitude, latitude = v12.equatorial_to_galactic(ra.ravel(), dec.ravel())
        membership = template_membership(longitude, latitude).reshape(
            n, n, template_count
        )
        total = membership.sum(axis=(0, 1), dtype=np.int64).astype(float)
        row_sum = membership.sum(axis=1, dtype=np.int64).astype(float)
        column_sum = membership.sum(axis=0, dtype=np.int64).astype(float)
        flat = membership.reshape(n * n, template_count).astype(np.float64)
        overlap = flat.T @ flat

        mean += total / n
        if n > 1:
            residual_cross = (
                overlap
                - (row_sum.T @ row_sum) / n
                - (column_sum.T @ column_sum) / n
                + np.outer(total, total) / (n * n)
            )
            covariance += residual_cross / (n - 1)

        if block_index % 100 == 0 or block_index + 1 == len(state.blocks):
            print(f"blocks {block_index + 1}/{len(state.blocks)}", flush=True)

    return mean, covariance


def gaussian_max_probability(
    correlation: np.ndarray,
    threshold: float,
    trials: int = 5_000_000,
    seed: int = 20260820,
) -> tuple[float, int, float]:
    symmetric = (correlation + correlation.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    factor = eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))
    rng = np.random.default_rng(seed)
    exceedances = 0
    completed = 0
    while completed < trials:
        count = min(100_000, trials - completed)
        draws = rng.standard_normal((count, correlation.shape[0])) @ factor.T
        exceedances += int(np.count_nonzero(draws.max(axis=1) >= threshold))
        completed += count
    return (exceedances + 1) / (trials + 1), exceedances, float(eigenvalues.min())


def main() -> None:
    observed_hash = sha256(EVENTS_PATH)
    if observed_hash != EXPECTED_EVENTS_SHA256:
        raise RuntimeError(f"V12 event input identity mismatch: {observed_hash}")

    v12 = load_v12()
    events = pd.read_csv(EVENTS_PATH)
    observed_membership = template_membership(
        events["gal_l_deg"].to_numpy(float), events["gal_b_deg"].to_numpy(float)
    )
    observed = observed_membership.sum(axis=0, dtype=np.int64)
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
    mean, covariance = exact_joint_moments(v12, state)
    standard_deviation = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    z = np.divide(
        observed - mean,
        standard_deviation,
        out=np.full_like(mean, np.nan),
        where=standard_deviation > 0,
    )
    correlation = covariance / np.outer(standard_deviation, standard_deviation)
    correlation = np.nan_to_num(correlation)
    np.fill_diagonal(correlation, 1.0)
    threshold = float(np.nanmax(z))
    family_p, exceedances, minimum_eigenvalue = gaussian_max_probability(
        correlation, threshold
    )

    table = pd.DataFrame(
        {
            "template": TEMPLATE_LABELS,
            "N": observed,
            "B_exact": mean,
            "sB_exact": standard_deviation,
            "Z_exact": z,
            "normal_one_sided_p": norm.sf(z),
        }
    ).sort_values("Z_exact", ascending=False)
    summary = {
        "template_definitions_verified_before_run": True,
        "template_count": len(TEMPLATE_LABELS),
        "observed_maximum_z": threshold,
        "leading_template": str(table.iloc[0]["template"]),
        "gaussian_joint_family_p": family_p,
        "gaussian_trials": 5_000_000,
        "gaussian_exceedances": exceedances,
        "correlation_minimum_eigenvalue_before_clipping": minimum_eigenvalue,
        "direct_randomization_validation_required_by_locked_rule": bool(family_p < 0.05),
        "events_sha256": observed_hash,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT_DIR / "uhe_galactic_plane_template_statistics.csv", index=False)
    np.savez_compressed(
        OUTPUT_DIR / "uhe_galactic_plane_template_moments.npz",
        template_labels=np.asarray(TEMPLATE_LABELS),
        observed=observed,
        mean=mean,
        covariance=covariance,
        correlation=correlation,
        z=z,
    )
    (OUTPUT_DIR / "uhe_galactic_plane_template_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(table.to_string(index=False))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
