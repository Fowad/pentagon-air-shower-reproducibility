#!/usr/bin/env python3
"""Build the compact numerical and identity ledger for the V12A handoff."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audit"
RESULTS = AUDIT / "results"
EXPECTED_EVENT_SHA = "f7f370ad96d69ef9a793472490b4f10f4deff4e6c62cd2211527f1c65e245a6b"
EXPECTED_V12_SHA = "853c690442f3f956661bb41707192cbb1a9438a6f7be17e13f85c48bb1ea3435"


def read_json(path: Path):
    return json.loads(path.read_text())


def checkpoint_identity(path: Path):
    with np.load(path, allow_pickle=False) as data:
        completed = int(data["completed"])
        if "identity_json" in data.files:
            identity = json.loads(str(data["identity_json"]))
            target = int(identity["trials"])
            event_sha = identity["event_sha256"]
            v12_sha = identity["v12_code_sha256"]
        else:
            target = int(data["target_trials"])
            event_sha = str(data["event_sha256"])
            v12_sha = str(data["v12_code_sha256"])
        value_shape = list(data["values"].shape)
    return {
        "file": path.name,
        "completed": completed,
        "target_trials": target,
        "value_shape": value_shape,
        "event_sha256": event_sha,
        "v12_code_sha256": v12_sha,
        "identity_pass": (
            completed == target == 100_000
            and event_sha == EXPECTED_EVENT_SHA
            and v12_sha == EXPECTED_V12_SHA
        ),
    }


def main() -> None:
    manifest = read_json(AUDIT / "v12" / "V12_RUN_MANIFEST.json")
    primary = read_json(AUDIT / "v12" / "results" / "primary_results_V12.json")
    catalogs = read_json(RESULTS / "catalog_cache_integrity_v12a.json")
    lhaaso = read_json(RESULTS / "lhaaso_37_stack_summary_v12a.json")
    p300 = read_json(RESULTS / "p300_exact_family_summary_v12a.json")
    lsi_phase = read_json(RESULTS / "lsi_phase_bins_100k_summary.json")
    lsi_interval = read_json(RESULTS / "lsi_interval1_100k_summary.json")
    lsi_grid = read_json(RESULTS / "lsi_interval_phase_grid_100k_summary.json")
    maxi = read_json(RESULTS / "maxi_j1820_interval2_100k_summary.json")

    checkpoints = [
        checkpoint_identity(RESULTS / "lsi_phase_bins_100k_checkpoint.npz"),
        checkpoint_identity(RESULTS / "lsi_interval1_100k_checkpoint.npz"),
        checkpoint_identity(RESULTS / "lsi_interval_phase_grid_100k_checkpoint.npz"),
        checkpoint_identity(RESULTS / "maxi_j1820_interval2_100k_checkpoint.npz"),
    ]
    integrity_pass = (
        manifest["stage"] == "complete"
        and manifest["code_sha256"] == EXPECTED_V12_SHA
        and manifest["accepted_events"] == 532_527
        and manifest["event_sha256"] == EXPECTED_EVENT_SHA
        and manifest["reference_realizations"] == 100_000
        and manifest["look_elsewhere_realizations"] == 100_000
        and primary["event_sha256"] == EXPECTED_EVENT_SHA
        and catalogs["all_downloaded_hashes_match"]
        and catalogs["required_failure_count"] == 0
        and all(item["identity_pass"] for item in checkpoints)
    )

    result = {
        "release_id": "PENTAGON-PEV-FOLLOWUP-V12A-2026-08-20",
        "integrity_pass": integrity_pass,
        "v12": {
            "release": manifest["release"],
            "code_sha256": manifest["code_sha256"],
            "accepted_events": manifest["accepted_events"],
            "event_sha256": manifest["event_sha256"],
            "run_count": primary["run_count"],
            "livetime_hours": primary["livetime_hours"],
            "reference_realizations": primary["reference_realizations"],
            "calibration_realizations": primary["calibration_realizations"],
            "blind_full_map_p": primary["primary_map_p"],
            "blind_six_interval_family_p": primary["six_interval_family_p"],
            "gamma_all_valid_full_p": primary["gamma_catalog_all_valid_full_p"],
            "gamma_all_valid_six_interval_family_p": primary[
                "gamma_catalog_all_valid_six_interval_family_p"
            ],
            "gamma_high_occupancy_full_p": primary[
                "gamma_catalog_high_occupancy_full_p"
            ],
            "gamma_high_occupancy_six_interval_family_p": primary[
                "gamma_catalog_high_occupancy_six_interval_family_p"
            ],
        },
        "catalog_cache": {
            key: catalogs[key]
            for key in (
                "catalog_family_count",
                "downloaded_snapshot_count",
                "embedded_published_table_count",
                "catalog_row_total",
                "required_failure_count",
                "all_downloaded_hashes_match",
            )
        },
        "lhaaso_uhe_37": {
            key: lhaaso[key]
            for key in (
                "full_archive_max_source",
                "full_archive_max_z",
                "full_archive_37_source_family_p",
                "equal_source_stouffer_z",
                "equal_source_stouffer_one_sided_p",
                "flux100_weighted_stouffer_z",
                "flux100_weighted_stouffer_one_sided_p",
                "family_gaussian_trials",
            )
        },
        "energy_matched_17": {
            key: p300[key]
            for key in (
                "full_observed_max_source",
                "full_observed_max_z",
                "full_gaussian_family_p",
                "interval_observed_max_source",
                "interval_observed_max_interval",
                "interval_observed_max_z",
                "interval_gaussian_family_p",
                "full_equal_source_stack_z",
                "interval_equal_source_six_interval_family_p",
                "interval_aggregate_count_six_interval_family_p",
                "gaussian_trials",
            )
        },
        "lsi_direct": {
            "standalone_interval_1_local_p": lsi_interval["one_sided_p"],
            "joint_grid_interval_1_local_p": lsi_grid["interval_tests"][0][
                "one_sided_p"
            ],
            "source_specific_six_interval_family_p": lsi_grid[
                "six_interval_family_p"
            ],
            "full_archive_local_p": lsi_phase["full_archive"]["one_sided_p"],
            "ten_phase_family_p": lsi_phase["ten_phase_bin_family_p"],
            "two_published_phase_window_family_p": lsi_phase[
                "two_published_window_family_p"
            ],
            "twelve_interval_published_window_family_p": lsi_grid[
                "twelve_interval_published_window_family_p"
            ],
            "sixty_interval_phase_cell_family_p": lsi_grid[
                "sixty_interval_phase_cell_family_p"
            ],
        },
        "maxi_j1820_interval_2": {
            "Z": maxi["Z"],
            "local_p": maxi["one_sided_p"],
            "parent_17_by_6_family_p": p300["interval_gaussian_family_p"],
        },
        "direct_100k_checkpoints": checkpoints,
        "bounded_conclusion": (
            "No statistically significant source-associated directional excess "
            "is resolved in the accepted five-detector event stream under the "
            "audited blind, catalog-directed, UHE/PeV, and externally motivated "
            "temporal tests. This is not a primary-gamma flux or source-inactivity "
            "exclusion."
        ),
    }
    output = RESULTS / "V12A_NUMERICAL_CROSSCHECK.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"integrity_pass": integrity_pass, "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
