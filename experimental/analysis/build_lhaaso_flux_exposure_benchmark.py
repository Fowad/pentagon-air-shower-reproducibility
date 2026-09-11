#!/usr/bin/env python3
"""Fold first-LHAASO catalog power laws through actual directional live time.

This is a flux/exposure benchmark only.  It does not apply a detector effective
area and therefore does not predict accepted events.  The catalog power law is
used exactly as published for the KM2A component and is explicitly treated as
an optimistic extrapolation where spectral curvature is unknown.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


BANDS_TEV = ((25.0, 100.0), (100.0, 300.0), (300.0, 1000.0), (1000.0, 1600.0))
REFERENCE_AREAS_M2 = (100.0, 1000.0, 10_000.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def power_law_integral(
    norm_50tev_1e16: np.ndarray,
    index: np.ndarray,
    low_tev: float,
    high_tev: float,
) -> np.ndarray:
    """Integral of N0*(E/50 TeV)^(-index), in photons cm^-2 s^-1."""

    norm = np.asarray(norm_50tev_1e16, dtype=float) * 1e-16
    gamma = np.asarray(index, dtype=float)
    if np.any(np.isclose(gamma, 1.0)):
        raise ValueError("index=1 requires the logarithmic special case")
    return (
        norm
        * 50.0
        / (gamma - 1.0)
        * (
            (low_tev / 50.0) ** (1.0 - gamma)
            - (high_tev / 50.0) ** (1.0 - gamma)
        )
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).parents[1]
    parser.add_argument(
        "--sources", type=Path, default=root / "audit" / "lhaaso_37_sources.csv"
    )
    parser.add_argument(
        "--exposure-summary",
        type=Path,
        default=Path(__file__).with_name("source_exposure")
        / "lhaaso_37_actual_live_exposure_summary.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).with_name("source_exposure"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.sources.resolve()
    exposure_path = args.exposure_summary.resolve()
    sources = pd.read_csv(source_path)
    exposure = pd.read_csv(exposure_path)[
        ["source", "visible_live_hours_theta_lt_45"]
    ]
    table = sources.merge(
        exposure, on="source", how="inner", validate="one_to_one"
    )
    visible_seconds = table.visible_live_hours_theta_lt_45.to_numpy(float) * 3600.0
    for low, high in BANDS_TEV:
        tag = f"{low:g}_{high:g}tev"
        integral = power_law_integral(
            table.flux_norm_50tev.to_numpy(float),
            table.spectral_index.to_numpy(float),
            low,
            high,
        )
        table[f"integral_flux_{tag}_cm-2_s-1"] = integral
        table[f"visible_fluence_{tag}_cm-2"] = integral * visible_seconds
        for area in REFERENCE_AREAS_M2:
            table[f"incident_photons_if_constant_{area:g}m2_{tag}"] = (
                integral * visible_seconds * area * 1e4
            )

    combined_integral = power_law_integral(
        table.flux_norm_50tev.to_numpy(float),
        table.spectral_index.to_numpy(float),
        25.0,
        1600.0,
    )
    table["integral_flux_25_1600tev_cm-2_s-1"] = combined_integral
    table["visible_fluence_25_1600tev_cm-2"] = combined_integral * visible_seconds
    for area in REFERENCE_AREAS_M2:
        table[f"incident_photons_if_constant_{area:g}m2_25_1600tev"] = (
            combined_integral * visible_seconds * area * 1e4
        )

    table = table.sort_values(
        "visible_fluence_300_1000tev_cm-2", ascending=False
    )
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "lhaaso_37_flux_live_time_benchmark.csv"
    json_path = output_dir / "LHAASO_FLUX_LIVE_TIME_BENCHMARK.json"
    table.to_csv(csv_path, index=False)
    leading = table.iloc[0]
    summary = {
        "protocol": "PENTAGON-LHAASO-FLUX-LIVE-TIME-BENCHMARK-V1",
        "source_spectrum": (
            "first-LHAASO KM2A catalog power law, dN/dE=N0*(E/50 TeV)^(-Gamma)"
        ),
        "interpretation_constraint": (
            "incident-photon benchmarks are not accepted-event predictions; "
            "the catalog warns that about one third of KM2A components are curved"
        ),
        "bands_tev": [list(band) for band in BANDS_TEV],
        "reference_areas_m2": list(REFERENCE_AREAS_M2),
        "source_count": int(len(table)),
        "leading_300_1000tev_source": str(leading.source),
        "leading_300_1000tev_visible_fluence_cm-2": float(
            leading["visible_fluence_300_1000tev_cm-2"]
        ),
        "leading_incident_photons_at_10000m2": float(
            leading["incident_photons_if_constant_10000m2_300_1000tev"]
        ),
        "sources": str(source_path),
        "sources_sha256": sha256(source_path),
        "exposure_summary": str(exposure_path),
        "exposure_summary_sha256": sha256(exposure_path),
        "output_csv": str(csv_path),
        "output_csv_sha256": sha256(csv_path),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    columns = [
        "source",
        "visible_live_hours_theta_lt_45",
        "visible_fluence_300_1000tev_cm-2",
        "incident_photons_if_constant_10000m2_300_1000tev",
    ]
    print(table[columns].head(12).to_string(index=False))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
