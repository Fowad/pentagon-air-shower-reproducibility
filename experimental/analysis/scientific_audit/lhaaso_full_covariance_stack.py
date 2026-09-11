#!/usr/bin/env python3
"""Exact full-archive covariance and stacking tests for 1LHAASO UHE sources."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from exact_permutation_moments import ROOT, load_v12, observed_counts, source_vectors


def exact_joint_moments(v12, state, sources):
    vectors = source_vectors(sources)
    count = len(sources)
    separations = np.degrees(
        np.arccos(np.clip(vectors @ vectors.T, -1.0, 1.0))
    )
    overlap_pairs = [
        (a, b)
        for a in range(count)
        for b in range(a, count)
        if a == b or separations[a, b] < 16.0 + 1.0e-12
    ]
    mean = np.zeros(count, dtype=float)
    covariance = np.zeros((count, count), dtype=float)
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
        ra = np.radians(ra.ravel())
        dec = np.radians(dec.ravel())
        sky_vectors = np.column_stack(
            (np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))
        )
        membership = ((sky_vectors @ vectors.T) >= cos_radius).reshape(n, n, count)
        total = membership.sum(axis=(0, 1), dtype=np.int64).astype(float)
        row_sum = membership.sum(axis=1, dtype=np.int64).astype(float)
        col_sum = membership.sum(axis=0, dtype=np.int64).astype(float)
        overlap = np.zeros((count, count), dtype=float)
        for a, b in overlap_pairs:
            value = float(np.count_nonzero(membership[:, :, a] & membership[:, :, b]))
            overlap[a, b] = value
            overlap[b, a] = value

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


def stack_result(observed, mean, covariance, source_weights):
    std = np.sqrt(np.diag(covariance))
    z = (observed - mean) / std
    correlation = covariance / np.outer(std, std)
    weights = np.asarray(source_weights, dtype=float)
    statistic = float(weights @ z / np.sqrt(weights @ correlation @ weights))
    return statistic, float(norm.sf(statistic))


def gaussian_max_probability(correlation, threshold, trials=5_000_000, seed=370001):
    eigenvalues, eigenvectors = np.linalg.eigh((correlation + correlation.T) / 2)
    factor = eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))
    rng = np.random.default_rng(seed)
    exceedances = 0
    for start in range(0, trials, 50_000):
        count = min(50_000, trials - start)
        values = rng.standard_normal((count, correlation.shape[0])) @ factor.T
        exceedances += int(np.count_nonzero(values.max(axis=1) >= threshold))
    return (
        float((1 + exceedances) / (trials + 1)),
        exceedances,
        float(eigenvalues.min()),
    )


def main():
    v12 = load_v12()
    sources = pd.read_csv(ROOT / "audit" / "lhaaso_37_sources.csv")
    events = pd.read_csv(ROOT / "audit" / "v12" / "reconstructed_events.csv.gz")
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
    observed = observed_counts(v12, events, sources)
    mean, covariance = exact_joint_moments(v12, state, sources)
    equal_z, equal_p = stack_result(observed, mean, covariance, np.ones(len(sources)))
    # The catalog's published 50-TeV differential-flux normalizations provide
    # an externally defined UHE-brightness proxy.  This is a sensitivity check,
    # not a detector-response-optimal weight because this array lacks gamma A_eff(E).
    flux_z, flux_p = stack_result(
        observed, mean, covariance, sources["flux_norm_50tev"].to_numpy(float)
    )
    # Convert the catalog power-law normalizations to 100 TeV.  This is the
    # preferred UHE brightness proxy for this follow-up and avoids weighting
    # the stack by a lower-energy normalization when the scientific question
    # concerns >100-TeV sources.
    flux100 = sources["flux_norm_50tev"].to_numpy(float) * np.power(
        2.0, -sources["spectral_index"].to_numpy(float)
    )
    flux100_z, flux100_p = stack_result(
        observed, mean, covariance, flux100
    )
    std = np.sqrt(np.diag(covariance))
    source_z = (observed - mean) / std
    correlation = covariance / np.outer(std, std)
    max_z = float(np.max(source_z))
    max_index = int(np.argmax(source_z))
    family_p, family_exceedances, correlation_min_eigenvalue = (
        gaussian_max_probability(correlation, max_z)
    )
    result = {
        "source_count": len(sources),
        "equal_source_stouffer_z": equal_z,
        "equal_source_stouffer_one_sided_p": equal_p,
        "flux50_weighted_stouffer_z": flux_z,
        "flux50_weighted_stouffer_one_sided_p": flux_p,
        "flux100_weighted_stouffer_z": flux100_z,
        "flux100_weighted_stouffer_one_sided_p": flux100_p,
        "full_archive_max_source": str(sources.iloc[max_index]["source"]),
        "full_archive_max_z": max_z,
        "full_archive_37_source_family_p": family_p,
        "full_archive_37_source_family_exceedances": family_exceedances,
        "family_gaussian_trials": 5_000_000,
        "correlation_min_eigenvalue": correlation_min_eigenvalue,
        "covariance_min_eigenvalue": float(np.linalg.eigvalsh(covariance).min()),
        "note": (
            "Exact permutation mean/covariance; Gaussian tail from 1279 independent "
            "four-hour blocks. The 100-TeV flux weighting is the preferred catalog-"
            "brightness sensitivity check. Neither flux weighting is detector-response "
            "optimal because gamma effective area versus energy is unavailable."
        ),
    }
    output_dir = ROOT / "audit" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "lhaaso_37_exact_full_covariance.npz",
        source=np.asarray(sources["source"], dtype=str),
        observed=observed,
        mean=mean,
        covariance=covariance,
    )
    (output_dir / "lhaaso_37_stack_summary_v12a.json").write_text(
        json.dumps(result, indent=2)
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
