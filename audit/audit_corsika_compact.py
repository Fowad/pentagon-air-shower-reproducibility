#!/usr/bin/env python3
"""Reaggregate and validate the compact CORSIKA R3 campaign products."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
COMPACT = ROOT / "corsika_compact" / "OVERNIGHT_STAGE2_R3_COMPACT_RESULT"
SUPPORT = ROOT / "simulation_support" / "PENTAGON_CORSIKA_PAPER_OVERNIGHT_STAGE2_R3"
OUTPUT = ROOT / "fresh_evidence" / "corsika_compact_reaggregation.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    support_manifest = SUPPORT / "PACKAGE_SHA256SUMS.txt"
    support_checked = 0
    support_failures = []
    for line in support_manifest.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        wanted, relative = line.split("  ", 1)
        path = SUPPORT / relative
        if not path.is_file():
            support_failures.append(f"missing: {relative}")
        elif sha256(path) != wanted:
            support_failures.append(f"hash mismatch: {relative}")
        support_checked += 1

    state = json.loads((COMPACT / "OVERNIGHT_STAGE2_STATE_R3.json").read_text())
    tasks = state["tasks"]
    manifests = sorted((COMPACT / "manifests").glob("*_manifest.json"))
    pilot_manifest = json.loads(
        (COMPACT / "ORIGINAL_STAGE2A_STAGE2A_PILOT_MANIFEST.json").read_text()
    )
    all_records = pd.read_csv(COMPACT / "OVERNIGHT_RESPONSE_ALL_RECORDS.csv")
    published_summary = pd.read_csv(COMPACT / "OVERNIGHT_RESPONSE_SUMMARY.csv")
    particle_audit = pd.read_csv(COMPACT / "PARTICLE_CLASS_DETECTOR_RESPONSE_AUDIT.csv")

    manifest_checks = []
    binary_hashes = set()
    source_hashes = set()
    for path in manifests:
        item = json.loads(path.read_text())
        label = str(item["label"])
        task = tasks[label]
        binary_hashes.add(item["binary_sha256"])
        source_hashes.add(item["minimal_source_sha256"])
        manifest_checks.append(
            {
                "label": label,
                "completed": task["status"] == "completed",
                "state_manifest_seed_match": int(item["seed"]) == int(task["seed"]),
                "state_manifest_ground_rows_match": int(item["ground_rows"])
                == int(task["ground_rows"]),
                "unit_particle_weights": bool(item["all_particle_weights_one"]),
            }
        )

    log_checks = []
    for label, task in tasks.items():
        log = COMPACT / "logs" / f"{label}_seed_{task['seed']}.log"
        text = log.read_text(errors="replace")
        primary_pattern = (
            rf"Primary: gamma, E={int(task['energy_gev'])} GeV, "
            rf"zenith={int(task['zenith_deg'])} deg, azimuth={int(task['azimuth_deg'])} deg"
        )
        log_checks.append(
            {
                "label": label,
                "exists": log.is_file(),
                "corsika_version_0p1p99": "This is CORSIKA8 0.1.99" in text,
                "declared_models": all(
                    token in text
                    for token in (
                        "SIBYLL-2.3d HE",
                        "UrQMD LE",
                        "PROPOSAL EM",
                        "SOPHIA resonance photonuclear",
                        "USStdBK atmosphere",
                    )
                ),
                "igrf13_site": (
                    "IGRF-13 field: year=2017 latitude=35.704 longitude=51.351 altitude=1200 m"
                    in text
                ),
                "primary_matches_state": re.search(primary_pattern, text) is not None,
                "energy_accounting_present": "Energy accounting (GeV):" in text,
            }
        )

    keys = ["energy_gev", "zenith_deg", "azimuth_deg", "pulse_window_ns", "threshold_mip"]
    regenerated_summary = (
        all_records.groupby(keys, as_index=False)
        .agg(
            n_fold_records=("effective_area_m2", "size"),
            n_simulation_labels=("simulation_label", "nunique"),
            accepted_core_points_min=("accepted_core_points", "min"),
            accepted_core_points_median=("accepted_core_points", "median"),
            accepted_core_points_max=("accepted_core_points", "max"),
            area_proxy_m2_min=("effective_area_m2", "min"),
            area_proxy_m2_median=("effective_area_m2", "median"),
            area_proxy_m2_max=("effective_area_m2", "max"),
            p68_deg_median=("p68_angular_error_deg", "median"),
        )
        .sort_values(keys)
        .reset_index(drop=True)
    )
    regenerated_summary.insert(
        1, "energy_pev", regenerated_summary["energy_gev"].to_numpy(float) / 1.0e6
    )
    regenerated_summary = regenerated_summary[published_summary.columns]
    expected = published_summary.sort_values(keys).reset_index(drop=True)
    columns_equal = list(regenerated_summary.columns) == list(expected.columns)
    summary_equal = columns_equal and all(
        np.allclose(
            regenerated_summary[column].to_numpy(float),
            expected[column].to_numpy(float),
            rtol=1e-12,
            atol=1e-12,
            equal_nan=True,
        )
        for column in expected.columns
        if column not in keys
    ) and all(
        np.allclose(
            regenerated_summary[column].to_numpy(float),
            expected[column].to_numpy(float),
            rtol=0,
            atol=0,
        )
        for column in keys
    )

    nominal = all_records[
        (all_records.zenith_deg == 20)
        & (all_records.azimuth_deg == 0)
        & (all_records.pulse_window_ns == 10)
        & (all_records.threshold_mip == 0.5)
    ]
    per_shower = (
        nominal.groupby(["energy_gev", "simulation_label"], as_index=False)
        .agg(
            area_response_m2=("effective_area_m2", "median"),
            area_8deg_m2=("source_aperture_effective_area_m2", "median"),
            r68_deg=("p68_angular_error_deg", "median"),
        )
        .sort_values(["energy_gev", "simulation_label"])
    )
    table_rows = []
    for energy, group in per_shower.groupby("energy_gev"):
        table_rows.append(
            {
                "energy_tev": float(energy / 1000.0),
                "showers": int(len(group)),
                "simulation_labels": group.simulation_label.tolist(),
                "area_response_m2_min": float(group.area_response_m2.min()),
                "area_response_m2_max": float(group.area_response_m2.max()),
                "area_8deg_m2_min": float(group.area_8deg_m2.min()),
                "area_8deg_m2_max": float(group.area_8deg_m2.max()),
                "r68_deg_min": float(group.r68_deg.min()),
                "r68_deg_max": float(group.r68_deg.max()),
            }
        )

    em = particle_audit.assign(
        electromagnetic_fraction=np.where(
            particle_audit.particle_class.isin(["e+", "e-", "gamma"]),
            particle_audit.fraction_of_sampled_detector_deposit,
            0.0,
        )
    )
    em = (
        em.groupby(["simulation_label", "energy_gev", "zenith_deg", "response_seed"])
        .electromagnetic_fraction.sum()
        .reset_index()
    )
    populated = em[(em.zenith_deg == 20) & (em.energy_gev >= 200_000)]

    z40 = all_records[all_records.zenith_deg == 40]
    z0 = all_records[all_records.zenith_deg == 0]
    z0_nominal = z0[(z0.pulse_window_ns == 10) & (z0.threshold_mip == 0.5)]

    missing_raw_parquet = not any(COMPACT.glob("**/particles.parquet"))
    raw_inventory_text = (COMPACT / "RAW_OUTPUTS_RETAIN_LOCALLY.txt").read_text()
    raw_inventory_count = raw_inventory_text.count("particles.parquet SHA256:")

    compact_hashes = {
        str(path.relative_to(COMPACT)): sha256(path)
        for path in sorted(COMPACT.glob("*.csv")) + sorted(COMPACT.glob("*.json"))
    }
    report = {
        "protocol": "PENTAGON-CORSIKA-COMPACT-REAGGREGATION-V1",
        "status": "archived_campaign_verified_fresh_shower_rerun_not_possible_from_compact_bundle",
        "campaign": {
            "r3_completed_tasks": len(tasks),
            "prior_pilot_showers": 1,
            "total_independent_corsika_showers_represented": len(tasks) + 1,
            "r3_elapsed_hours": float(state["elapsed_hours"]),
            "all_tasks_completed": all(item["status"] == "completed" for item in tasks.values()),
            "manifests": len(manifests),
            "manifest_checks_all_pass": all(
                all(value for key, value in row.items() if key != "label") for row in manifest_checks
            ),
            "log_checks_all_pass": all(
                all(value for key, value in row.items() if key != "label") for row in log_checks
            ),
            "corsika_version": "0.1.99",
            "locked_commit_claimed_by_package": "601fe3036725af876d6ab46a7cfe15f0c0987be5",
            "binary_sha256_values": sorted(binary_hashes),
            "minimal_source_sha256_values": sorted(source_hashes),
            "pilot_unit_particle_weights": bool(pilot_manifest["all_particle_weights_one"]),
            "manifest_checks": manifest_checks,
            "log_checks": log_checks,
        },
        "support_package_integrity": {
            "manifest_entries_checked": support_checked,
            "failures": support_failures,
            "passes": not support_failures,
        },
        "reaggregation": {
            "all_record_rows": int(len(all_records)),
            "expected_rows_from_19_showers_x_3_replicas_x_3_windows_x_4_thresholds": 684,
            "summary_rows": int(len(published_summary)),
            "regenerated_summary_exact_within_1e_minus_12": bool(summary_equal),
            "nominal_table_20deg_0p5mip_10ns": table_rows,
            "all_40deg_accepted_core_points_zero": bool((z40.accepted_core_points == 0).all()),
            "all_40deg_fivefold_core_points_zero": bool((z40.fivefold_core_points == 0).all()),
            "z0_nominal_nonzero_at_100_200_300tev": bool(
                (z0_nominal.groupby("energy_gev").accepted_core_points.max() > 0).all()
            ),
            "electromagnetic_deposit_fraction": {
                "all_replica_min": float(em.electromagnetic_fraction.min()),
                "all_replica_median": float(em.electromagnetic_fraction.median()),
                "all_replica_max": float(em.electromagnetic_fraction.max()),
                "fraction_of_all_replicas_above_0p98": float(
                    (em.electromagnetic_fraction > 0.98).mean()
                ),
                "z_le_20_min": float(
                    em[em.zenith_deg <= 20].electromagnetic_fraction.min()
                ),
                "populated_200_300tev_z20_min": float(
                    populated.electromagnetic_fraction.min()
                ),
                "populated_200_300tev_z20_median": float(
                    populated.electromagnetic_fraction.median()
                ),
            },
        },
        "raw_shower_reproducibility_limit": {
            "raw_particles_parquet_present": not missing_raw_parquet,
            "raw_inventory_entries": raw_inventory_count,
            "corsika_binary_present": any(ROOT.glob("**/c8_air_shower")),
            "corsika_source_tree_present": any(ROOT.glob("**/corsika8-source")),
            "meaning": (
                "The compact archive supports hash/configuration/log validation and exact reaggregation, "
                "but excludes the raw ground-particle Parquet files and the executable/source tree needed "
                "for a new first-principles shower generation or independent refold."
            ),
        },
        "top_level_compact_hashes": compact_hashes,
        "scientific_interpretation": (
            "The archived compact products exactly reproduce the manuscript's response table and qualitative "
            "energy/zenith claims. They are a sparse diagnostic campaign, not a precision effective-area "
            "calibration; fresh CORSIKA shower regeneration remains unverified in this environment."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
