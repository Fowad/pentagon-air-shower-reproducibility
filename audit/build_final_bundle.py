#!/usr/bin/env python3
"""Assemble the compact, independently verifiable V14 audit deliverable."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DELIVERABLES = ROOT / "deliverables"
BUNDLE = DELIVERABLES / "Pentagon_Array_V14_Final_Reproducibility_Audit"
ZIP_BASE = DELIVERABLES / "Pentagon_Array_V14_Final_Reproducibility_Audit_Bundle"
INTACT_V12 = Path("/workspace/scratch/c675a63d2681/audit_work/v12_inner")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_file(source: Path, relative: str) -> None:
    destination = BUNDLE / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def main() -> None:
    if BUNDLE.exists() or ZIP_BASE.with_suffix(".zip").exists():
        raise SystemExit("Refusing to overwrite an existing final bundle")
    BUNDLE.mkdir(parents=True)

    copy_file(ROOT / "V14_FINAL_REPRODUCIBILITY_AUDIT.md", "V14_FINAL_REPRODUCIBILITY_AUDIT.md")
    copy_file(
        ROOT / "V14_DISCREPANCY_AND_LIMITATION_REGISTER.csv",
        "V14_DISCREPANCY_AND_LIMITATION_REGISTER.csv",
    )
    shutil.copytree(ROOT / "fresh_evidence", BUNDLE / "fresh_evidence")

    script_names = [
        "audit_corsika_compact.py",
        "audit_events_coordinates_exposure.py",
        "build_final_bundle.py",
        "build_v14_claim_ledger.py",
        "compare_geant4_response.py",
        "run_geant4_parallel_audit.py",
        "timing_jitter_sensitivity.py",
        "verify_source_families.py",
    ]
    for name in script_names:
        copy_file(ROOT / "tools" / name, f"audit_scripts/{name}")

    for name in (
        "mc_checkpoint_reference_V12.npz",
        "mc_checkpoint_calibration_V12.npz",
        "V12_RUN_MANIFEST.json",
        "run_summary.csv",
    ):
        copy_file(INTACT_V12 / name, f"authoritative_checkpoints/v12/{name}")
    copy_file(
        INTACT_V12 / "results" / "primary_results_V12.json",
        "authoritative_checkpoints/v12/primary_results_V12.json",
    )

    results = ROOT / "recovered" / "scientific_audit_outputs" / "results"
    directed_names = [
        "V12A_NUMERICAL_CROSSCHECK.json",
        "lhaaso_37_exact_full_covariance.npz",
        "lhaaso_37_stack_summary_v12a.json",
        "p300_exact_family_covariance.npz",
        "p300_exact_family_summary_v12a.json",
        "lsi_interval1_100k_checkpoint.npz",
        "lsi_interval1_100k_summary.json",
        "lsi_phase_bins_100k_checkpoint.npz",
        "lsi_phase_bins_100k_summary.json",
        "lsi_interval_phase_grid_100k_checkpoint.npz",
        "lsi_interval_phase_grid_100k_summary.json",
        "lsi_published_windows_100k_checkpoint.npz",
        "lsi_orbital_phase_summary.json",
        "maxi_j1820_interval2_100k_checkpoint.npz",
        "maxi_j1820_interval2_100k_summary.json",
    ]
    for name in directed_names:
        copy_file(results / name, f"authoritative_checkpoints/directed_sources/{name}")

    shutil.copytree(
        ROOT / "recovered" / "aperture_robustness",
        BUNDLE / "authoritative_checkpoints" / "aperture_robustness",
    )
    shutil.copytree(
        ROOT / "recovered" / "source_exposure",
        BUNDLE / "authoritative_checkpoints" / "source_exposure",
    )
    exposure_expected = (
        ROOT
        / "recovered"
        / "code_and_reference_package"
        / "experimental"
        / "expected"
        / "source_exposure"
    )
    copy_file(
        exposure_expected / "lhaaso_37_flux_live_time_benchmark.csv",
        "authoritative_checkpoints/source_exposure/lhaaso_37_flux_live_time_benchmark.csv",
    )
    copy_file(
        exposure_expected / "LHAASO_FLUX_LIVE_TIME_BENCHMARK.json",
        "authoritative_checkpoints/source_exposure/LHAASO_FLUX_LIVE_TIME_BENCHMARK.json",
    )

    shutil.copytree(
        ROOT / "geant4_fresh" / "geant4_response",
        BUNDLE / "fresh_geant4" / "geant4_response",
    )
    shutil.copytree(
        ROOT / "corsika_compact" / "OVERNIGHT_STAGE2_R3_COMPACT_RESULT",
        BUNDLE / "archived_corsika_compact" / "OVERNIGHT_STAGE2_R3_COMPACT_RESULT",
    )
    copy_file(
        ROOT
        / "simulation_support"
        / "PENTAGON_CORSIKA_PAPER_OVERNIGHT_STAGE2_R3"
        / "PACKAGE_SHA256SUMS.txt",
        "archived_corsika_compact/SUPPORT_PACKAGE_SHA256SUMS.txt",
    )
    copy_file(
        ROOT / "recovered" / "code_and_reference_package" / "EXPECTED_IDENTITIES.json",
        "package_metadata/EXPECTED_IDENTITIES.json",
    )
    copy_file(
        ROOT / "recovered" / "code_and_reference_package" / "REPRODUCIBILITY_SPEC.md",
        "package_metadata/REPRODUCIBILITY_SPEC.md",
    )
    copy_file(
        ROOT / "prior_v13_evidence" / "independent_crosschecks.json",
        "prior_independent_evidence/independent_crosschecks.json",
    )
    copy_file(
        ROOT / "prior_v13_evidence" / "supporting_crosschecks.json",
        "prior_independent_evidence/supporting_crosschecks.json",
    )

    files = sorted(path for path in BUNDLE.rglob("*") if path.is_file())
    manifest = {
        "protocol": "PENTAGON-V14-FINAL-AUDIT-BUNDLE-V1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "file_count_before_manifests": len(files),
        "notes": [
            "The raw DAQ ZIP, accepted-event cache, manuscript PDF, and original code ZIP are not duplicated.",
            "Their identities and the audit scope are recorded in V14_FINAL_REPRODUCIBILITY_AUDIT.md.",
            "The fresh GEANT4 event basis and all numerical checkpoints directly used in this audit are included.",
        ],
    }
    (BUNDLE / "BUNDLE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    files = sorted(path for path in BUNDLE.rglob("*") if path.is_file())
    hash_lines = [
        f"{sha256(path)}  {path.relative_to(BUNDLE).as_posix()}" for path in files
    ]
    (BUNDLE / "SHA256SUMS.txt").write_text("\n".join(hash_lines) + "\n")

    zip_path = Path(
        shutil.make_archive(
            str(ZIP_BASE),
            "zip",
            root_dir=DELIVERABLES,
            base_dir=BUNDLE.name,
        )
    )
    result = {
        "bundle_directory": str(BUNDLE),
        "bundle_zip": str(zip_path),
        "zip_bytes": zip_path.stat().st_size,
        "zip_sha256": sha256(zip_path),
        "files_in_bundle": sum(1 for path in BUNDLE.rglob("*") if path.is_file()),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
