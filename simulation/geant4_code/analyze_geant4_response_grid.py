#!/usr/bin/env python3
"""Validate and summarize the GEANT4 particle-response grid."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MIP_MPV_MEV = 3.57852
THRESHOLDS_MIP = (0.25, 0.5, 1.0, 2.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = successes / trials
    denom = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denom
    half = z * np.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / denom
    return float(center - half), float(center + half)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--grid-dir",
        type=Path,
        default=Path(__file__).with_name("geant4_response") / "particle_grid_v1",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grid_dir = args.grid_dir.resolve()
    manifest_path = grid_dir / "GEANT4_PARTICLE_RESPONSE_GRID_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    records = manifest["records"]
    if len(records) != 270 or any(row["status"] == "failed" for row in records):
        raise RuntimeError("response-grid manifest is incomplete")

    wide_rows: list[dict] = []
    long_rows: list[dict] = []
    for record in records:
        path = Path(record["output"])
        if not path.is_file():
            path = grid_dir / path.name
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise RuntimeError(f"response-grid hash mismatch: {path}")
        with np.load(path, allow_pickle=False) as data:
            deposits = np.asarray(data["deposited_energy_mev"], dtype=float)
            particle = str(data["particle"])
            energy = float(data["kinetic_energy_mev"])
            angle = float(data["angle_deg"])
            include_steel = bool(data["include_steel"])
            geant4_version = int(data["geant4_version_number"])
        if not include_steel:
            raise RuntimeError(f"main response grid unexpectedly lacks steel: {path}")
        quantile_probabilities = (0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
        quantiles = np.quantile(deposits, quantile_probabilities)
        wide = {
            "particle": particle,
            "kinetic_energy_mev": energy,
            "angle_deg": angle,
            "events": len(deposits),
            "include_steel": include_steel,
            "geant4_version_number": geant4_version,
            "nonzero_probability": float(np.mean(deposits > 0)),
            "mean_deposit_mev": float(np.mean(deposits)),
            "std_deposit_mev": float(np.std(deposits, ddof=1)),
        }
        for probability, value in zip(quantile_probabilities, quantiles):
            wide[f"q{probability:g}_deposit_mev"] = float(value)
        for threshold in THRESHOLDS_MIP:
            cutoff = threshold * MIP_MPV_MEV
            successes = int(np.count_nonzero(deposits >= cutoff))
            low, high = wilson_interval(successes, len(deposits))
            probability = successes / len(deposits)
            wide[f"p_ge_{threshold:g}mip"] = probability
            long_rows.append(
                {
                    "particle": particle,
                    "kinetic_energy_mev": energy,
                    "angle_deg": angle,
                    "threshold_mip": threshold,
                    "threshold_mev": cutoff,
                    "events": len(deposits),
                    "successes": successes,
                    "response_probability": probability,
                    "wilson95_low": low,
                    "wilson95_high": high,
                    "include_steel": include_steel,
                }
            )
        wide_rows.append(wide)

    wide_table = pd.DataFrame(wide_rows).sort_values(
        ["particle", "kinetic_energy_mev", "angle_deg"]
    )
    long_table = pd.DataFrame(long_rows).sort_values(
        ["particle", "threshold_mip", "angle_deg", "kinetic_energy_mev"]
    )
    wide_path = grid_dir / "geant4_particle_response_summary.csv"
    long_path = grid_dir / "geant4_particle_threshold_probabilities.csv"
    wide_table.to_csv(wide_path, index=False)
    long_table.to_csv(long_path, index=False)

    colors = {0.25: "#2166ac", 0.5: "#67a9cf", 1.0: "#ef8a62", 2.0: "#b2182b"}
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9), constrained_layout=True)
    for axis, particle, title in zip(
        axes[:2], ("gamma", "e-"), ("Ground photons", "Ground electrons")
    ):
        subset = long_table[(long_table.particle == particle) & (long_table.angle_deg == 20)]
        for threshold in THRESHOLDS_MIP:
            curve = subset[subset.threshold_mip == threshold]
            axis.plot(
                curve.kinetic_energy_mev,
                curve.response_probability,
                marker="o",
                ms=3.5,
                lw=1.4,
                color=colors[threshold],
                label=f"{threshold:g} MIP",
            )
        axis.set_xscale("log")
        axis.set_ylim(-0.03, 1.03)
        axis.set_xlabel("Kinetic energy (MeV)")
        axis.set_ylabel("Single-particle station-response probability")
        axis.set_title(f"{title}, 20° incidence")
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)

    representative = {
        "gamma": (100.0, "100 MeV photon", "#762a83"),
        "e-": (100.0, "100 MeV electron", "#1b7837"),
        "mu-": (10_000.0, "10 GeV muon", "#d95f02"),
    }
    for particle, (energy, label, color) in representative.items():
        curve = long_table[
            (long_table.particle == particle)
            & (long_table.kinetic_energy_mev == energy)
            & (long_table.threshold_mip == 1.0)
        ]
        axes[2].plot(
            curve.angle_deg,
            curve.response_probability,
            marker="o",
            lw=1.5,
            color=color,
            label=label,
        )
    axes[2].set_ylim(-0.03, 1.03)
    axes[2].set_xlabel("Incidence angle (deg)")
    axes[2].set_ylabel("Response probability at 1 MIP")
    axes[2].set_title("Path-length dependence")
    axes[2].grid(alpha=0.25)
    axes[2].legend(frameon=False, fontsize=8)
    figure_path = grid_dir / "geant4_particle_response_validation.png"
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)

    summary = {
        "status": "complete",
        "protocol": manifest["protocol"],
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "files": len(records),
        "total_events": int(wide_table.events.sum()),
        "mip_mpv_mev": MIP_MPV_MEV,
        "thresholds_mip": list(THRESHOLDS_MIP),
        "geant4_version_numbers": sorted(
            int(value) for value in wide_table.geant4_version_number.unique()
        ),
        "wide_csv": str(wide_path),
        "wide_csv_sha256": sha256(wide_path),
        "long_csv": str(long_path),
        "long_csv_sha256": sha256(long_path),
        "figure": str(figure_path),
        "figure_sha256": sha256(figure_path),
    }
    summary_path = grid_dir / "geant4_particle_response_audit.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
