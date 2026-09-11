#!/usr/bin/env python3
"""Exact LS I +61 303 source-circle tests in externally published phase bins."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from exact_permutation_moments import ROOT, load_v12


SOURCE_RA_DEG = 40.1319
SOURCE_DEC_DEG = 61.2293
PHASE_EPOCH_MJD = 43366.275
ORBIT_DAYS = 26.496


def phase_masks(events):
    times = pd.to_datetime(events["utc_time"], utc=True, format="mixed")
    mjd = times.astype("int64").to_numpy(float) / 1.0e9 / 86400.0 + 40587.0
    phase = np.mod(mjd - PHASE_EPOCH_MJD, ORBIT_DAYS) / ORBIT_DAYS
    labels = ["full", "published_25_100TeV_phase_0.3_0.6", "published_gt100TeV_phase_0.6_0.2"]
    masks = [
        np.ones(len(events), dtype=bool),
        (phase >= 0.3) & (phase < 0.6),
        (phase >= 0.6) | (phase < 0.2),
    ]
    for index in range(10):
        low = index / 10.0
        high = (index + 1) / 10.0
        labels.append(f"phase_{low:.1f}_{high:.1f}")
        masks.append((phase >= low) & (phase < high))
    return labels, np.asarray(masks), phase


def exact_joint_moments(v12, state, labels, masks):
    q = len(labels)
    mean = np.zeros(q, dtype=float)
    covariance = np.zeros((q, q), dtype=float)
    ra0 = math.radians(SOURCE_RA_DEG)
    dec0 = math.radians(SOURCE_DEC_DEG)
    source_vector = np.array(
        [math.cos(dec0) * math.cos(ra0), math.cos(dec0) * math.sin(ra0), math.sin(dec0)]
    )
    cos_radius = math.cos(math.radians(8.0)) - 1.0e-15

    for block_index, indices in enumerate(state.blocks):
        indices = np.asarray(indices, dtype=np.int64)
        n = len(indices)
        hour_angle = state.hour_angle_deg[indices][:, None]
        dec_mean = np.broadcast_to(
            state.declination_mean_deg[indices][:, None], (n, n)
        )
        lst = np.broadcast_to(state.observed_lst_deg[indices][None, :], (n, n))
        angles = [
            np.broadcast_to(angle[indices][None, :], (n, n))
            for angle in state.precession
        ]
        ra, dec = v12.mean_of_date_to_icrs_with_angles(
            (lst - hour_angle) % 360.0, dec_mean, *angles
        )
        ra = np.radians(ra)
        dec = np.radians(dec)
        membership = (
            np.cos(dec) * np.cos(ra) * source_vector[0]
            + np.cos(dec) * np.sin(ra) * source_vector[1]
            + np.sin(dec) * source_vector[2]
        ) >= cos_radius
        statistics = membership[:, :, None] & masks[:, indices].T[None, :, :]
        total = statistics.sum(axis=(0, 1), dtype=np.int64).astype(float)
        row_sum = statistics.sum(axis=1, dtype=np.int64).astype(float)
        col_sum = statistics.sum(axis=0, dtype=np.int64).astype(float)
        flat = statistics.reshape(n * n, q).astype(float)
        overlap = flat.T @ flat
        mean += total / n
        if n > 1:
            covariance += (
                overlap
                - (row_sum.T @ row_sum) / n
                - (col_sum.T @ col_sum) / n
                + np.outer(total, total) / (n * n)
            ) / (n - 1)
        if block_index % 100 == 0 or block_index + 1 == len(state.blocks):
            print(f"blocks {block_index + 1}/{len(state.blocks)}", flush=True)
    return mean, covariance


def gaussian_max_family(correlation, threshold, trials, seed):
    eigenvalues, eigenvectors = np.linalg.eigh((correlation + correlation.T) / 2)
    factor = eigenvectors @ np.diag(np.sqrt(np.clip(eigenvalues, 0.0, None)))
    rng = np.random.default_rng(seed)
    exceedances = 0
    for start in range(0, trials, 100_000):
        count = min(100_000, trials - start)
        draws = rng.standard_normal((count, len(correlation))) @ factor.T
        exceedances += int(np.count_nonzero(draws.max(axis=1) >= threshold))
    return (1 + exceedances) / (trials + 1), exceedances, float(eigenvalues.min())


def main():
    v12 = load_v12()
    events = pd.read_csv(ROOT / "audit" / "v12" / "reconstructed_events.csv.gz")
    labels, masks, phase = phase_masks(events)
    inside = v12._exact_circle_membership(
        events["ra_deg"].to_numpy(float),
        events["dec_deg"].to_numpy(float),
        SOURCE_RA_DEG,
        SOURCE_DEC_DEG,
        8.0,
    )
    observed = np.asarray([np.count_nonzero(inside & mask) for mask in masks])
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
    mean, covariance = exact_joint_moments(v12, state, labels, masks)
    std = np.sqrt(np.diag(covariance))
    z = (observed - mean) / std
    p = norm.sf(z)
    correlation = covariance / np.outer(std, std)

    two = np.array([1, 2])
    two_threshold = float(z[two].max())
    two_p, two_exceed, two_min_eig = gaussian_max_family(
        correlation[np.ix_(two, two)], two_threshold, 5_000_000, 161303
    )
    ten = np.arange(3, 13)
    ten_threshold = float(z[ten].max())
    ten_p, ten_exceed, ten_min_eig = gaussian_max_family(
        correlation[np.ix_(ten, ten)], ten_threshold, 5_000_000, 26496
    )

    table = pd.DataFrame(
        {
            "test": labels,
            "N": observed,
            "B_exact": mean,
            "sB_exact": std,
            "Z_exact": z,
            "normal_one_sided_p": p,
        }
    )
    summary = {
        "source": "LS I +61 303",
        "ra_deg": SOURCE_RA_DEG,
        "dec_deg": SOURCE_DEC_DEG,
        "phase_epoch_mjd": PHASE_EPOCH_MJD,
        "period_days": ORBIT_DAYS,
        "published_phase_window_family_max_z": two_threshold,
        "published_phase_window_family_p": two_p,
        "published_phase_window_family_exceedances": two_exceed,
        "ten_phase_bin_max_z": ten_threshold,
        "ten_phase_bin_family_p": ten_p,
        "ten_phase_bin_family_exceedances": ten_exceed,
        "published_window_correlation_min_eigenvalue": two_min_eig,
        "ten_bin_correlation_min_eigenvalue": ten_min_eig,
        "gaussian_trials": 5_000_000,
        "note": (
            "Exact V12 permutation means/covariances. Gaussian family tails use the "
            "joint limit from 1279 independent four-hour blocks. The two broad phase "
            "windows were fixed by the LHAASO energy-resolved result before examining "
            "this array."
        ),
    }
    output_dir = ROOT / "audit" / "results"
    table.to_csv(output_dir / "lsi_orbital_phase_exact.csv", index=False)
    (output_dir / "lsi_orbital_phase_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(table.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
