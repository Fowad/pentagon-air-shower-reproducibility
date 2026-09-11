#!/usr/bin/env python3
"""Regenerate the supplied GEANT4 response basis with process-level parallelism.

The authoritative generator deliberately executes each cell in a fresh serial
Geant4 process.  This audit driver preserves that behavior, but schedules the
independent cells concurrently.  Once all cells are fresh and validated, it
invokes the authoritative aggregator, which revalidates every NPZ and writes
the canonical manifests and CSV summaries.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def load_generator(path: Path):
    spec = importlib.util.spec_from_file_location("pentagon_geant4_generator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import generator: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    generator_path = args.generator.resolve()
    workspace = args.workspace.resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    g = load_generator(generator_path)

    response_root = workspace / "geant4_response"
    mip_dir = response_root / "mip_calibration_v1"
    grid_dir = response_root / "particle_grid_v1"
    high_dir = response_root / "high_energy_em_extension_v1"
    seeds = g.archived_seed_lookup()

    tasks: list[dict] = [
        {
            "kind": "mip",
            "output": mip_dir / "vertical_muminus_10GeV_20000.npz",
            "log": mip_dir / "vertical_muminus_10GeV_20000.log",
            "particle": "mu-",
            "energy": 10000.0,
            "angle": 0.0,
            "seed": g.MIP_SEED,
            "events": g.MIP_EVENTS,
            "include_steel": True,
        }
    ]

    grid_cells = [
        (particle, energy, angle)
        for particle in g.EM_PARTICLES
        for energy in g.EM_ENERGIES_MEV
        for angle in g.ANGLES
    ] + [
        (particle, energy, angle)
        for particle in g.OTHER_PARTICLES
        for energy in g.OTHER_ENERGIES_MEV
        for angle in g.ANGLES
    ]
    for particle, energy, angle in grid_cells:
        events = g.EVENTS_EM if particle in g.EM_PARTICLES else g.EVENTS_OTHER
        seed = seeds.get(
            (particle, float(energy), float(angle), False),
            g.stable_seed(particle, energy, angle, "steel"),
        )
        stem = f"{g.particle_tag(particle)}_E{g.energy_tag(energy)}MeV_A{angle:g}deg"
        tasks.append(
            {
                "kind": "grid",
                "output": grid_dir / f"{stem}.npz",
                "log": grid_dir / "logs" / f"{stem}.log",
                "particle": particle,
                "energy": energy,
                "angle": angle,
                "seed": seed,
                "events": events,
                "include_steel": True,
            }
        )

    high_cells = [
        (particle, energy, angle, False)
        for particle in g.EM_PARTICLES
        for energy in g.HIGH_EM_ENERGIES_MEV
        for angle in g.ANGLES
    ] + [
        ("gamma", energy, angle, True)
        for energy in g.HIGH_EM_ENERGIES_MEV
        for angle in g.ANGLES
    ]
    for particle, energy, angle, bare in high_cells:
        seed = seeds.get(
            (particle, float(energy), float(angle), bare),
            g.stable_seed(particle, energy, angle, "bare" if bare else "steel"),
        )
        prefix = "bare" if bare else "steel"
        stem = f"{prefix}_{g.particle_tag(particle)}_E{g.energy_tag(energy)}MeV_A{angle:g}deg"
        tasks.append(
            {
                "kind": "high",
                "output": high_dir / f"{stem}.npz",
                "log": high_dir / "logs" / f"{stem}.log",
                "particle": particle,
                "energy": energy,
                "angle": angle,
                "seed": seed,
                "events": g.EVENTS_HIGH_EM,
                "include_steel": not bare,
            }
        )

    started = time.time()
    records: list[dict] = []

    def execute(task: dict) -> dict:
        result = g.run_cell_process(
            task["output"],
            task["log"],
            task["particle"],
            task["energy"],
            task["angle"],
            task["seed"],
            task["events"],
            task["include_steel"],
        )
        return {
            **task,
            "output": str(task["output"]),
            "log": str(task["log"]),
            **result,
            "sha256": g.sha256(task["output"]),
        }

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(execute, task) for task in tasks]
        for index, future in enumerate(as_completed(futures), 1):
            record = future.result()
            records.append(record)
            print(
                f"[{index:03d}/{len(tasks)}] {record['kind']} "
                f"{record['particle']} {record['energy']:g}MeV "
                f"{record['angle']:g}deg",
                flush=True,
            )

    generated_manifest = {
        "protocol": "PENTAGON-GEANT4-PARALLEL-INDEPENDENT-RERUN-V1",
        "authoritative_generator": str(generator_path),
        "authoritative_generator_sha256": g.sha256(generator_path),
        "workers": args.workers,
        "tasks": len(tasks),
        "fresh_tasks": sum(not record["reused"] for record in records),
        "reused_tasks": sum(record["reused"] for record in records),
        "runtime_seconds": time.time() - started,
        "records": sorted(records, key=lambda row: row["output"]),
    }
    manifest_path = response_root / "PARALLEL_RERUN_MANIFEST.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(generated_manifest, indent=2, sort_keys=True) + "\n")

    subprocess.run(
        [sys.executable, str(generator_path), "--all", "--workspace", str(workspace)],
        check=True,
    )
    print(json.dumps({k: generated_manifest[k] for k in (
        "workers", "tasks", "fresh_tasks", "reused_tasks", "runtime_seconds"
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
