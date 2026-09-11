#!/usr/bin/env python3
"""Fast static QA for the public repository; no raw data or heavy simulation needed."""
from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
PROHIBITED_SUFFIXES = {".parquet", ".root", ".h5", ".hdf5"}
PROHIBITED_NAMES = {".DS_Store"}


def fail(message: str) -> None:
    raise SystemExit("QA FAILED: " + message)


def main() -> None:
    py_files = sorted(ROOT.rglob("*.py"))
    for path in py_files:
        try:
            compile(path.read_text(encoding="utf-8"), str(path), "exec")
        except Exception as exc:
            fail(f"Python parse error in {path.relative_to(ROOT)}: {exc}")

    for path in ROOT.rglob("*.json"):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            fail(f"JSON parse error in {path.relative_to(ROOT)}: {exc}")

    for path in [ROOT / "environment.yml", ROOT / "CITATION.cff"]:
        try:
            yaml.safe_load(path.read_text(encoding="utf-8"))
        except Exception as exc:
            fail(f"YAML parse error in {path.relative_to(ROOT)}: {exc}")

    for path in ROOT.rglob("*.csv"):
        try:
            with path.open(newline="", encoding="utf-8", errors="replace") as handle:
                next(csv.reader(handle), None)
        except Exception as exc:
            fail(f"CSV read error in {path.relative_to(ROOT)}: {exc}")

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.name in PROHIBITED_NAMES:
            fail(f"OS metadata present: {path.relative_to(ROOT)}")
        if path.suffix.lower() in PROHIBITED_SUFFIXES:
            fail(f"Raw/bulky data type present: {path.relative_to(ROOT)}")
        if path.suffix.lower() in {".npz", ".npy"}:
            fail(f"Generated checkpoint array present in public repo: {path.relative_to(ROOT)}")

    bash = shutil.which("bash")
    commands = sorted(ROOT.rglob("*.command"))
    if bash:
        for path in commands:
            proc = subprocess.run([bash, "-n", str(path)], capture_output=True, text=True)
            if proc.returncode:
                fail(f"Shell syntax error in {path.relative_to(ROOT)}: {proc.stderr}")

    print(f"QA OK: {len(py_files)} Python files, {len(commands)} shell launchers; no raw/bulky data found.")


if __name__ == "__main__":
    main()
