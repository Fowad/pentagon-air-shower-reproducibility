#!/usr/bin/env python3
"""Check the GEANT4 electromagnetic-response plateau at 30 and 100 GeV."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np


ENERGIES_MEV = (30_000.0, 100_000.0)
ANGLES_DEG = (0, 20, 40, 60, 70)
PARTICLES = ("gamma", "e-", "e+")


@dataclass(frozen=True)
class Task:
    particle: str
    energy_mev: float
    angle_deg: float
    bare: bool
    events: int
    seed: int
    output: str
    log: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_seed(particle: str, energy: float, angle: float, bare: bool) -> int:
    label = f"PENTAGON-G4-HIGH-EM-V1|{particle}|{energy:.12g}|{angle:.12g}|bare={bare}"
    return 1 + int.from_bytes(hashlib.sha256(label.encode()).digest()[:4], "big") % 2_000_000_000


def validate(task: Task) -> bool:
    path = Path(task.output)
    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as data:
            return (
                str(data["particle"]) == task.particle
                and float(data["kinetic_energy_mev"]) == task.energy_mev
                and float(data["angle_deg"]) == task.angle_deg
                and bool(data["include_steel"]) == (not task.bare)
                and int(data["seed"]) == task.seed
                and len(data["deposited_energy_mev"]) == task.events
            )
    except Exception:
        return False


def run_one(task: Task, driver: Path) -> dict:
    if validate(task):
        return {**asdict(task), "status": "reused", "sha256": sha256(Path(task.output))}
    command = [
        sys.executable,
        str(driver),
        "--particle",
        task.particle,
        "--energy-mev",
        str(task.energy_mev),
        "--angle-deg",
        str(task.angle_deg),
        "--events",
        str(task.events),
        "--seed",
        str(task.seed),
        "--output",
        task.output,
    ]
    if task.bare:
        command.append("--bare")
    Path(task.output).parent.mkdir(parents=True, exist_ok=True)
    Path(task.log).parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with Path(task.log).open("wb") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    runtime = time.monotonic() - started
    if result.returncode or not validate(task):
        return {
            **asdict(task),
            "status": "failed",
            "returncode": result.returncode,
            "runtime_seconds": runtime,
        }
    return {
        **asdict(task),
        "status": "completed",
        "returncode": result.returncode,
        "runtime_seconds": runtime,
        "sha256": sha256(Path(task.output)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--driver",
        type=Path,
        default=Path(__file__).with_name("pentagon_scintillator_geant4.py"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).with_name("geant4_response")
        / "high_energy_em_extension_v1",
    )
    parser.add_argument("--events", type=int, default=2_000)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.events <= 0 or args.workers <= 0:
        raise ValueError("events and workers must be positive")
    driver = args.driver.resolve()
    output_dir = args.output_dir.resolve()
    tasks = []
    for bare, particles in ((False, PARTICLES), (True, ("gamma",))):
        geometry = "bare" if bare else "steel"
        for particle in particles:
            particle_tag = particle.replace("-", "minus").replace("+", "plus")
            for energy in ENERGIES_MEV:
                for angle in ANGLES_DEG:
                    stem = f"{geometry}_{particle_tag}_E{energy:g}MeV_A{angle:g}deg".replace(".", "p")
                    tasks.append(
                        Task(
                            particle=particle,
                            energy_mev=energy,
                            angle_deg=angle,
                            bare=bare,
                            events=args.events,
                            seed=deterministic_seed(particle, energy, angle, bare),
                            output=str(output_dir / f"{stem}.npz"),
                            log=str(output_dir / "logs" / f"{stem}.log"),
                        )
                    )
    records = []
    started = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_one, task, driver): task for task in tasks}
        for number, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            records.append(row)
            print(
                f"[{number:02d}/{len(tasks):02d}] {row['status']:9s} "
                f"{row['particle']} {row['energy_mev']/1000:g} GeV "
                f"{row['angle_deg']:g} deg bare={row['bare']}",
                flush=True,
            )
    records.sort(key=lambda r: (r["bare"], r["particle"], r["energy_mev"], r["angle_deg"]))
    manifest = {
        "protocol": "PENTAGON-GEANT4-HIGH-ENERGY-EM-EXTENSION-V1",
        "driver": str(driver),
        "driver_sha256": sha256(driver),
        "energies_mev": list(ENERGIES_MEV),
        "angles_deg": list(ANGLES_DEG),
        "events_per_cell": args.events,
        "started_unix": started,
        "finished_unix": time.time(),
        "records": records,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "GEANT4_HIGH_ENERGY_EM_EXTENSION_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    failures = [row for row in records if row["status"] == "failed"]
    print(f"manifest={manifest_path} failures={len(failures)}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
