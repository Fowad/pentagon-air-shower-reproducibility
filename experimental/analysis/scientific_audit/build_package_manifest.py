#!/usr/bin/env python3
"""Write a deterministic SHA-256 manifest for one deliverable directory."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / "PACKAGE_MANIFEST.json"
    files = []
    for path in sorted(root.rglob("*")):
        if (
            not path.is_file()
            or path == output
            or "__pycache__" in path.parts
        ):
            continue
        files.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest = {
        "release_id": "PENTAGON-PEV-FOLLOWUP-V12A-2026-08-20",
        "manifest_scope": "All package payload files except this manifest",
        "file_count": len(files),
        "payload_bytes": sum(item["bytes"] for item in files),
        "files": files,
    }
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({key: value for key, value in manifest.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
