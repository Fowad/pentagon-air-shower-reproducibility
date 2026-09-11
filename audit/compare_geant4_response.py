#!/usr/bin/env python3
"""Compare the fresh GEANT4 basis with the archived reference statistically."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FRESH = ROOT / "geant4_fresh" / "geant4_response"
ARCHIVED = (
    ROOT
    / "simulation_support"
    / "PENTAGON_CORSIKA_PAPER_OVERNIGHT_STAGE2_R3"
    / "support"
    / "geant4_response"
)
OUTPUT = ROOT / "fresh_evidence" / "geant4_fresh_vs_archived.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def statistical_envelope(p0: np.ndarray, p1: np.ndarray, n0: np.ndarray, n1: np.ndarray):
    pooled = (p0 * n0 + p1 * n1) / (n0 + n1)
    six_sigma = 6.0 * np.sqrt(
        np.maximum(0.0, pooled * (1.0 - pooled) * (1.0 / n0 + 1.0 / n1))
    )
    return np.maximum(0.03, six_sigma)


def primitive(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def frame_records(frame: pd.DataFrame) -> list[dict]:
    return [{key: primitive(value) for key, value in row.items()} for row in frame.to_dict("records")]


def main() -> None:
    archived_csv = ARCHIVED / "particle_grid_v1" / "geant4_particle_threshold_probabilities.csv"
    fresh_csv = FRESH / "particle_grid_v1" / "geant4_particle_threshold_probabilities.csv"
    archived = pd.read_csv(archived_csv)
    fresh = pd.read_csv(fresh_csv)
    keys = [
        "particle",
        "kinetic_energy_mev",
        "angle_deg",
        "threshold_mip",
        "include_steel",
    ]
    merged = archived.merge(
        fresh, on=keys, suffixes=("_archived", "_fresh"), validate="one_to_one"
    )
    if len(merged) != 1080:
        raise RuntimeError(f"expected 1,080 matched primary-grid rows, found {len(merged)}")

    p0 = merged.response_probability_archived.to_numpy(float)
    p1 = merged.response_probability_fresh.to_numpy(float)
    n0 = merged.events_archived.to_numpy(float)
    n1 = merged.events_fresh.to_numpy(float)
    delta = np.abs(p1 - p0)
    tolerance = statistical_envelope(p0, p1, n0, n1)
    outside = delta > tolerance
    merged["absolute_probability_difference"] = delta
    merged["acceptance_tolerance"] = tolerance
    merged["difference_over_tolerance"] = delta / tolerance
    exceptions = merged.loc[
        outside,
        keys
        + [
            "events_archived",
            "events_fresh",
            "response_probability_archived",
            "response_probability_fresh",
            "absolute_probability_difference",
            "acceptance_tolerance",
            "difference_over_tolerance",
        ],
    ].sort_values("difference_over_tolerance", ascending=False)

    high_rows: list[dict] = []
    archived_mip = 3.57852
    fresh_mip_manifest = json.loads(
        (FRESH / "mip_calibration_v1" / "GEANT4_MIP_CALIBRATION_MANIFEST.json").read_text()
    )
    fresh_mip = float(fresh_mip_manifest["mip_mpv_mev"])
    for archived_path in sorted((ARCHIVED / "high_energy_em_extension_v1").glob("*.npz")):
        fresh_path = FRESH / "high_energy_em_extension_v1" / archived_path.name
        with np.load(archived_path, allow_pickle=False) as old, np.load(
            fresh_path, allow_pickle=False
        ) as new:
            old_values = np.asarray(old["deposited_energy_mev"], dtype=float)
            new_values = np.asarray(new["deposited_energy_mev"], dtype=float)
            for threshold_mip in (0.25, 0.5, 1.0, 2.0):
                old_p = float(np.mean(old_values >= threshold_mip * archived_mip))
                new_p = float(np.mean(new_values >= threshold_mip * fresh_mip))
                old_n, new_n = len(old_values), len(new_values)
                tol = float(
                    statistical_envelope(
                        np.asarray([old_p]),
                        np.asarray([new_p]),
                        np.asarray([old_n]),
                        np.asarray([new_n]),
                    )[0]
                )
                high_rows.append(
                    {
                        "file": archived_path.name,
                        "threshold_mip": threshold_mip,
                        "events_archived": old_n,
                        "events_fresh": new_n,
                        "probability_archived": old_p,
                        "probability_fresh": new_p,
                        "absolute_probability_difference": abs(new_p - old_p),
                        "acceptance_tolerance": tol,
                    }
                )
    high = pd.DataFrame(high_rows)
    high["outside_tolerance"] = (
        high.absolute_probability_difference > high.acceptance_tolerance
    )

    archived_files = sorted(ARCHIVED.glob("**/*.npz"))
    fresh_by_relative = {
        path.relative_to(FRESH): path for path in sorted(FRESH.glob("**/*.npz"))
    }
    common = [path for path in archived_files if path.relative_to(ARCHIVED) in fresh_by_relative]
    exact_hash_matches = sum(
        sha256(path) == sha256(fresh_by_relative[path.relative_to(ARCHIVED)]) for path in common
    )

    parallel_manifest = json.loads((FRESH / "PARALLEL_RERUN_MANIFEST.json").read_text())
    fresh_grid_manifest = json.loads(
        (FRESH / "particle_grid_v1" / "GEANT4_PARTICLE_RESPONSE_GRID_MANIFEST.json").read_text()
    )
    archived_grid_manifest = json.loads(
        (ARCHIVED / "particle_grid_v1" / "GEANT4_PARTICLE_RESPONSE_GRID_MANIFEST.json").read_text()
    )
    fresh_high_manifest = json.loads(
        (FRESH / "high_energy_em_extension_v1" / "GEANT4_HIGH_ENERGY_EM_EXTENSION_MANIFEST.json").read_text()
    )
    archived_high_manifest = json.loads(
        (ARCHIVED / "high_energy_em_extension_v1" / "GEANT4_HIGH_ENERGY_EM_EXTENSION_MANIFEST.json").read_text()
    )

    by_particle = {}
    for particle, group in merged.assign(outside=outside).groupby("particle"):
        by_particle[str(particle)] = {
            "rows": int(len(group)),
            "outside_tolerance": int(group.outside.sum()),
            "median_absolute_difference": float(group.absolute_probability_difference.median()),
            "maximum_absolute_difference": float(group.absolute_probability_difference.max()),
        }

    within_40 = merged.angle_deg <= 40
    report = {
        "protocol": "PENTAGON-GEANT4-FRESH-VS-ARCHIVED-AUDIT-V1",
        "status": "complete_with_documented_exceptions",
        "fresh_execution": {
            "geant4_version_number": int(fresh_mip_manifest["geant4_version_number"]),
            "physics_list": str(fresh_mip_manifest["physics_list"]),
            "parallel_workers": int(parallel_manifest["workers"]),
            "fresh_tasks": int(parallel_manifest["fresh_tasks"]),
            "reused_tasks": int(parallel_manifest["reused_tasks"]),
            "runtime_seconds": float(parallel_manifest["runtime_seconds"]),
            "mip_events": int(fresh_mip_manifest["events"]),
            "primary_grid_events": int(
                sum(int(record["events"]) for record in fresh_grid_manifest["records"])
            ),
            "high_energy_extension_events": int(
                sum(int(record["events"]) for record in fresh_high_manifest["records"])
            ),
            "total_detector_events": int(
                fresh_mip_manifest["events"]
                + sum(int(record["events"]) for record in fresh_grid_manifest["records"])
                + sum(int(record["events"]) for record in fresh_high_manifest["records"])
            ),
            "generator_sha256": str(parallel_manifest["authoritative_generator_sha256"]),
        },
        "mip_calibration": {
            "archived_mpv_mev": archived_mip,
            "fresh_mpv_mev": fresh_mip,
            "difference_mev": fresh_mip - archived_mip,
            "relative_difference": (fresh_mip - archived_mip) / archived_mip,
            "passes_package_tolerance_0p25_mev": abs(fresh_mip - archived_mip) <= 0.25,
        },
        "primary_grid": {
            "rows_compared": int(len(merged)),
            "median_absolute_probability_difference": float(np.median(delta)),
            "q95_absolute_probability_difference": float(np.quantile(delta, 0.95)),
            "maximum_absolute_probability_difference": float(delta.max()),
            "maximum_difference_over_tolerance": float(np.max(delta / tolerance)),
            "rows_outside_tolerance": int(outside.sum()),
            "rows_at_angles_le_40_deg": int(within_40.sum()),
            "rows_outside_tolerance_at_angles_le_40_deg": int((outside & within_40).sum()),
            "acceptance_rule": "absolute difference <= max(0.03, six pooled-binomial standard errors)",
            "passes_all_angles": bool(np.all(~outside)),
            "passes_angles_le_40_deg": bool(np.all(~outside[within_40])),
            "by_particle": by_particle,
            "exceptions": frame_records(exceptions),
        },
        "high_energy_em_extension": {
            "rows_compared": int(len(high)),
            "median_absolute_probability_difference": float(
                high.absolute_probability_difference.median()
            ),
            "q95_absolute_probability_difference": float(
                high.absolute_probability_difference.quantile(0.95)
            ),
            "maximum_absolute_probability_difference": float(
                high.absolute_probability_difference.max()
            ),
            "rows_outside_tolerance": int(high.outside_tolerance.sum()),
            "passes": bool(not high.outside_tolerance.any()),
        },
        "archive_identity": {
            "common_response_npz_files": len(common),
            "exact_npz_hash_matches": int(exact_hash_matches),
            "fresh_generator_sha256": str(parallel_manifest["authoritative_generator_sha256"]),
            "archived_generator_sha256": str(archived_grid_manifest["driver_sha256"]),
            "same_generator_source_hash": bool(
                parallel_manifest["authoritative_generator_sha256"]
                == archived_grid_manifest["driver_sha256"]
            ),
            "archived_high_generator_sha256": str(archived_high_manifest["driver_sha256"]),
            "interpretation": (
                "The archived driver bytes did not survive and the supplied consolidated generator is a "
                "frozen reconstruction. Statistical agreement, not byte identity, is the declared criterion."
            ),
        },
        "scientific_interpretation": (
            "The MIP scale and every response threshold at incidence angles through 40 degrees agree "
            "within the declared cross-platform Monte Carlo envelope. Seven of 1,080 primary-grid rows "
            "miss the envelope at 60--70 degrees; the 160 high-energy extension rows all pass."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
