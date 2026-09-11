#!/usr/bin/env python3
"""Verify every downloaded V12 catalog snapshot against its run manifest."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "audit" / "catalog_cache_upload" / "catalog_cache_V12"
MANIFEST = CACHE / "catalog_snapshot_manifest_V12.csv"
V12_CODE = (
    ROOT
    / "audit"
    / "handoff"
    / "PENTAGON_V12_NEW_WORK_HANDOFF_PACKAGE"
    / "Pentagon_Array_Directional_Analysis_V12.py"
)
OUTPUT = ROOT / "audit" / "results" / "catalog_cache_integrity_v12a.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    table = pd.read_csv(MANIFEST)
    records = []
    for row in table.itertuples(index=False):
        snapshot = str(row.snapshot)
        record = {
            "catalog": str(row.catalog),
            "snapshot": snapshot,
            "rows": int(row.rows),
            "manifest_sha256": str(row.sha256),
        }
        if snapshot == "downloaded":
            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(row.catalog)).strip("_")
            candidates = sorted(CACHE.glob(f"{safe}.*"))
            if len(candidates) != 1:
                raise RuntimeError(
                    f"Expected one cached snapshot for {row.catalog!r}; found {candidates}"
                )
            actual = sha256(candidates[0])
            record.update(
                {
                    "file": candidates[0].name,
                    "actual_sha256": actual,
                    "hash_match": actual == str(row.sha256),
                }
            )
        else:
            record.update(
                {
                    "file": None,
                    "actual_sha256": None,
                    "hash_match": None,
                    "note": "Published table embedded in the definitive V12 code, not a downloaded cache file.",
                }
            )
        records.append(record)

    downloaded = [item for item in records if item["snapshot"] == "downloaded"]
    result = {
        "release_id": "PENTAGON-PEV-FOLLOWUP-V12A-2026-08-20",
        "catalog_manifest": MANIFEST.name,
        "catalog_manifest_sha256": sha256(MANIFEST),
        "v12_code_sha256": sha256(V12_CODE),
        "catalog_family_count": len(records),
        "downloaded_snapshot_count": len(downloaded),
        "embedded_published_table_count": len(records) - len(downloaded),
        "catalog_row_total": int(table["rows"].sum()),
        "required_failure_count": int(
            ((table["required"] == True) & (table["status"] != "loaded")).sum()
        ),
        "all_downloaded_hashes_match": all(item["hash_match"] for item in downloaded),
        "records": records,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "records"}, indent=2))


if __name__ == "__main__":
    main()
