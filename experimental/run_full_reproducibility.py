#!/usr/bin/env python3
"""Public-facing end-to-end experimental reproduction driver.

The 44 original DAQ text files are intentionally NOT distributed in this
repository. For byte-stable archival catalog provenance, pass the separately
retained frozen catalog cache with --catalog-cache-zip. Without it, the V12
analysis queries the official catalog services and records new snapshot hashes.

Default trial counts match the archival scientific run and are expensive.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE = HERE / "analysis"
FROZEN = HERE / "frozen_inputs"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def run(command: list[object], *, cwd: Path | None, log: Path) -> None:
    command = [str(x) for x in command]
    print("[RUN]", " ".join(command), flush=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("\n[RUN] " + " ".join(command) + "\n")
        handle.flush()
        proc = subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT)
    if proc.returncode:
        raise RuntimeError(f"Command failed with rc={proc.returncode}; see {log}")


def write_manifest(root: Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "FINAL_FILE_MANIFEST.json":
            rows.append({
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    (root / "FINAL_FILE_MANIFEST.json").write_text(
        json.dumps({"created_unix": time.time(), "files": rows}, indent=2) + "\n",
        encoding="utf-8",
    )


def stage_audit_layout(audit_root: Path, v12_out: Path) -> None:
    handoff = audit_root / "audit/handoff/PENTAGON_V12_NEW_WORK_HANDOFF_PACKAGE"
    v12_results = audit_root / "audit/v12/results"
    generic_results = audit_root / "audit/results"
    handoff.mkdir(parents=True, exist_ok=True)
    v12_results.mkdir(parents=True, exist_ok=True)
    generic_results.mkdir(parents=True, exist_ok=True)

    shutil.copy2(CODE / "Pentagon_Array_Directional_Analysis_V12.py",
                 handoff / "Pentagon_Array_Directional_Analysis_V12.py")
    shutil.copy2(FROZEN / "lhaaso_37_sources.csv", audit_root / "audit/lhaaso_37_sources.csv")
    shutil.copy2(FROZEN / "pevatron_sources_working.csv", audit_root / "audit/pevatron_sources_working.csv")
    shutil.copy2(v12_out / "reconstructed_events.csv.gz", audit_root / "audit/v12/reconstructed_events.csv.gz")
    shutil.copy2(v12_out / "run_summary.csv", audit_root / "audit/v12/run_summary.csv")

    package = v12_out / "Pentagon_Array_V12_RESULT_PACKAGE.zip"
    with zipfile.ZipFile(package) as archive:
        for name in archive.namelist():
            if name.startswith("results/") and not name.endswith("/"):
                destination = audit_root / "audit/v12" / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(archive.read(name))
    shutil.copytree(CODE / "scientific_audit", audit_root / "scientific_audit", dirs_exist_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Directory containing exactly the 44 original DAQ .txt files.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="New/empty output directory for the reproduction.")
    parser.add_argument("--catalog-cache-zip", type=Path, default=None,
                        help="Optional frozen catalog_cache_V12.zip from the archival package.")
    parser.add_argument("--reference-realizations", type=int, default=100000)
    parser.add_argument("--calibration-realizations", type=int, default=100000)
    parser.add_argument("--targeted-realizations", type=int, default=10000)
    parser.add_argument("--phase-trials", type=int, default=100000)
    parser.add_argument("--gaussian-family-trials", type=int, default=5000000)
    parser.add_argument("--aperture-trials", type=int, default=2000)
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not data_dir.is_dir():
        raise FileNotFoundError(data_dir)
    raw_files = sorted(data_dir.glob("*.txt"))
    if len(raw_files) != 44:
        raise RuntimeError(f"Expected exactly 44 DAQ .txt files; found {len(raw_files)}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(f"Output directory must be new or empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    v12_out = output_dir / "V12_MAIN"
    v12_out.mkdir()
    audit_root = output_dir / "downstream_reproducibility"
    log = output_dir / "MASTER_REPRODUCIBILITY_LOG.txt"

    if args.catalog_cache_zip is not None:
        cache_zip = args.catalog_cache_zip.expanduser().resolve()
        if not cache_zip.is_file():
            raise FileNotFoundError(cache_zip)
        with zipfile.ZipFile(cache_zip) as archive:
            archive.extractall(v12_out)
    else:
        print("[INFO] No frozen catalog cache supplied; official catalog services will be queried.")

    # Main 44-run reconstruction and definitive reference/calibration ensembles.
    run([
        sys.executable,
        CODE / "Pentagon_Array_Directional_Analysis_V12.py",
        "--data-dir", data_dir,
        "--output-dir", v12_out,
        "--reference-realizations", args.reference_realizations,
        "--look-elsewhere-realizations", args.calibration_realizations,
        "--targeted-source-realizations", args.targeted_realizations,
    ], cwd=HERE, log=log)

    # Downstream archival products used by the final manuscript audit.
    stage_audit_layout(audit_root, v12_out)
    sa = audit_root / "scientific_audit"
    run([sys.executable, sa / "lhaaso_full_covariance_stack.py"], cwd=audit_root, log=log)
    run([sys.executable, sa / "pevatron_family_covariance.py"], cwd=audit_root, log=log)
    run([sys.executable, sa / "validate_lsi_phase_bins_permutation.py",
         "--trials", args.phase_trials, "--workers", args.workers], cwd=audit_root, log=log)
    run([sys.executable, sa / "validate_lsi_interval_phase_grid_permutation.py",
         "--trials", args.phase_trials, "--workers", args.workers], cwd=audit_root, log=log)
    run([sys.executable, CODE / "archive_gaussian_family_nulls.py",
         "--root", audit_root, "--trials", args.gaussian_family_trials], cwd=HERE, log=log)
    run([sys.executable, CODE / "run_aperture_robustness_repro.py",
         "--root", audit_root, "--trials", args.aperture_trials,
         "--workers", args.workers], cwd=HERE, log=log)

    exposure = output_dir / "source_exposure"
    run([sys.executable, CODE / "build_source_exposure_kernel.py",
         "--run-summary", v12_out / "run_summary.csv",
         "--sources", FROZEN / "lhaaso_37_sources.csv",
         "--time-step-seconds", 30,
         "--zenith-bin-deg", 0.25,
         "--maximum-zenith-deg", 45,
         "--output-dir", exposure], cwd=HERE, log=log)
    run([sys.executable, CODE / "build_lhaaso_flux_exposure_benchmark.py",
         "--sources", FROZEN / "lhaaso_37_sources.csv",
         "--exposure-summary", exposure / "lhaaso_37_actual_live_exposure_summary.csv",
         "--output-dir", exposure], cwd=HERE, log=log)

    snapshot = output_dir / "EXECUTING_CODE_SNAPSHOT"
    shutil.copytree(HERE, snapshot, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.npz", "*.parquet"))
    write_manifest(output_dir)
    print("Complete:", output_dir)


if __name__ == "__main__":
    main()
