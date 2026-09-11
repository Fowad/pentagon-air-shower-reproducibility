#!/usr/bin/env python3
"""Isolate the directional smearing induced by the campaign's 1.85 ns jitter."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "fresh_evidence" / "electronics_jitter_sensitivity.json"
C_VACUUM_M_PER_S = 299_792_458.0
C_AIR_M_PER_S = 299_702_547.0
COORDINATES_M = np.asarray(
    [
        [-2.7586, -3.6756],
        [-4.3133, 1.0717],
        [-0.2780, 4.0196],
        [3.7707, 1.0942],
        [2.2376, -3.6617],
    ],
    dtype=float,
)


def simulate(theta_deg: float, normal_draws: np.ndarray, azimuth_rad: np.ndarray, sigma_ns: float):
    theta = np.deg2rad(theta_deg)
    truth_horizontal = np.sin(theta) * np.column_stack(
        [np.cos(azimuth_rad), np.sin(azimuth_rad)]
    )
    perfect_times = truth_horizontal @ COORDINATES_M.T / C_VACUUM_M_PER_S
    times = perfect_times + normal_draws * sigma_ns * 1e-9
    baselines = COORDINATES_M[1:] - COORDINATES_M[0]
    inverse = np.linalg.pinv(baselines)
    relative_distance = (times[:, 1:] - times[:, [0]]) * C_VACUUM_M_PER_S
    reconstructed_horizontal = relative_distance @ inverse.T
    horizontal_norm = np.linalg.norm(reconstructed_horizontal, axis=1)
    direction_valid = np.isfinite(horizontal_norm) & (horizontal_norm <= 1.0)
    reconstructed_theta_deg = np.full(len(times), np.nan)
    reconstructed_theta_deg[direction_valid] = np.rad2deg(
        np.arcsin(horizontal_norm[direction_valid])
    )
    reconstructed_vertical = np.sqrt(
        np.maximum(0.0, 1.0 - np.sum(reconstructed_horizontal**2, axis=1))
    )
    cosine = (
        np.sum(reconstructed_horizontal * truth_horizontal, axis=1)
        + reconstructed_vertical * np.cos(theta)
    )
    separation_deg = np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))
    accepted = direction_valid & (reconstructed_theta_deg < 45.0)
    accepted_errors = separation_deg[accepted]
    return {
        "theta_true_deg": theta_deg,
        "electronics_sigma_ns": sigma_ns,
        "trials": int(len(times)),
        "direction_valid_fraction": float(direction_valid.mean()),
        "theta_lt_45_acceptance_fraction": float(accepted.mean()),
        "median_angular_error_deg": float(np.median(accepted_errors)),
        "r68_angular_error_deg": float(np.quantile(accepted_errors, 0.68)),
        "r95_angular_error_deg": float(np.quantile(accepted_errors, 0.95)),
    }


def main() -> None:
    trial_count = 500_000
    sigmas_ns = (0.0, 0.5, 1.0, 1.85, 2.5)
    results = []
    for theta_deg in (0.0, 20.0, 40.0):
        rng = np.random.default_rng(20_260_902 + int(theta_deg))
        azimuth_rad = rng.uniform(0.0, 2.0 * np.pi, trial_count)
        normal_draws = rng.normal(0.0, 1.0, size=(trial_count, 5))
        for sigma_ns in sigmas_ns:
            results.append(simulate(theta_deg, normal_draws, azimuth_rad, sigma_ns))

    speed_effect = []
    ratio = C_VACUUM_M_PER_S / C_AIR_M_PER_S
    for theta_deg in (0.0, 20.0, 40.0):
        reconstructed = np.rad2deg(
            np.arcsin(np.clip(ratio * np.sin(np.deg2rad(theta_deg)), -1.0, 1.0))
        )
        speed_effect.append(
            {
                "theta_true_deg": theta_deg,
                "theta_reconstructed_deg": float(reconstructed),
                "deterministic_error_deg": float(reconstructed - theta_deg),
            }
        )

    nominal = [row for row in results if row["electronics_sigma_ns"] == 1.85]
    report = {
        "protocol": "PENTAGON-ELECTRONICS-JITTER-SENSITIVITY-V1",
        "method": (
            "Perfect planar fronts at random azimuths were sampled on the exact five-station geometry; "
            "independent Gaussian station-time errors were then added before the campaign's pseudoinverse "
            "direction reconstruction and theta<45 degree cut."
        ),
        "trial_count_per_theta_and_sigma": trial_count,
        "detector_coordinates_m": COORDINATES_M.tolist(),
        "reconstruction_speed_m_per_s": C_VACUUM_M_PER_S,
        "campaign_default_electronics_sigma_ns": 1.85,
        "results": results,
        "nominal_1p85ns_results": nominal,
        "vacuum_vs_air_speed_only": {
            "air_speed_m_per_s": C_AIR_M_PER_S,
            "vacuum_speed_m_per_s": C_VACUUM_M_PER_S,
            "results": speed_effect,
            "interpretation": "The speed choice contributes at most about 0.014 degree through 40 degrees and is negligible here.",
        },
        "scientific_interpretation": (
            "The campaign's 1.85 ns Gaussian station jitter alone produces roughly 8 degree r68 error "
            "for ideal planar fronts. It therefore supplies a large fraction of the reported near-11-degree "
            "simulation containment and must be disclosed and empirically justified before that simulation "
            "is called independent confirmation of the measured angular resolution."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
