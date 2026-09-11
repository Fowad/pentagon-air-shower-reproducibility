#!/usr/bin/env python3
"""Run and audit the predeclared GEANT4 ground-particle response grid."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


EM_ENERGIES_MEV = (0.5, 1, 3, 10, 30, 100, 300, 1_000, 3_000, 10_000)
OTHER_ENERGIES_MEV = (300, 1_000, 3_000, 10_000, 30_000, 100_000)
ANGLES_DEG = (0, 20, 40, 60, 70)
EM_PARTICLES = ("gamma", "e-", "e+")
OTHER_PARTICLES = ("mu-", "proton", "neutron", "pi+")


@dataclass(frozen=True)
class Task:
    particle: str
    energy_mev: float
    angle_deg: float
    events: int
    seed: int
    output: str
    log: str


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def energy_tag(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def make_seed(particle: str, energy: float, angle: float) -> int:
    payload = f"PENTAGON-G4-V1|{particle}|{energy:.12g}|{angle:.12g}".encode()
    return 1 + int.from_bytes(hashlib.sha256(payload).digest()[:4], "big") % 2_000_000_000


def build_tasks(output_dir: Path, events_em: int, events_other: int) -> list[Task]:
    tasks: list[Task] = []
    for particles, energies, events in (
        (EM_PARTICLES, EM_ENERGIES_MEV, events_em),
        (OTHER_PARTICLES, OTHER_ENERGIES_MEV, events_other),
    ):
        for particle in particles:
            particle_tag = particle.replace("+", "plus").replace("-", "minus")
            for energy in energies:
                for angle in ANGLES_DEG:
                    stem = f"{particle_tag}_E{energy_tag(energy)}MeV_A{angle:g}deg"
                    tasks.append(
                        Task(
                            particle=particle,
                            energy_mev=energy,
                            angle_deg=angle,
                            events=events,
                            seed=make_seed(particle, energy, angle),
                            output=str(output_dir / f"{stem}.npz"),
                            log=str(output_dir / "logs" / f"{stem}.log"),
                        )
                    )
    return tasks


def validate_existing(task: Task) -> bool:
    path = Path(task.output)
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        import numpy as np

        with np.load(path, allow_pickle=False) as data:
            return (
                len(data["deposited_energy_mev"]) == task.events
                and str(data["particle"]) == task.particle
                and float(data["kinetic_energy_mev"]) == task.energy_mev
                and float(data["angle_deg"]) == task.angle_deg
                and int(data["seed"]) == task.seed
            )
    except Exception:
        return False


def run_one(task: Task, driver: Path) -> dict:
    output = Path(task.output)
    log = Path(task.log)
    output.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    if validate_existing(task):
        return {
            **asdict(task),
            "status": "reused",
            "runtime_seconds": 0.0,
            "bytes": output.stat().st_size,
            "sha256": file_sha256(output),
        }
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
        str(output),
    ]
    start = time.monotonic()
    with log.open("wb") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    runtime = time.monotonic() - start
    if result.returncode != 0 or not validate_existing(task):
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
        "bytes": output.stat().st_size,
        "sha256": file_sha256(output),
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
        default=Path(__file__).with_name("geant4_response") / "particle_grid_v1",
    )
    parser.add_argument("--events-em", type=int, default=3_000)
    parser.add_argument("--events-other", type=int, default=2_000)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.events_em <= 0 or args.events_other <= 0 or args.workers <= 0:
        raise ValueError("event counts and worker count must be positive")
    driver = args.driver.resolve()
    if not driver.is_file():
        raise FileNotFoundError(driver)
    output_dir = args.output_dir.resolve()
    tasks = build_tasks(output_dir, args.events_em, args.events_other)
    print(f"tasks={len(tasks)} workers={args.workers} output={output_dir}")
    if args.dry_run:
        print(json.dumps([asdict(task) for task in tasks[:10]], indent=2))
        return

    started = time.time()
    records: list[dict] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_one, task, driver): task for task in tasks}
        for number, future in enumerate(as_completed(futures), start=1):
            record = future.result()
            records.append(record)
            print(
                f"[{number:03d}/{len(tasks):03d}] {record['status']:9s} "
                f"{record['particle']:7s} E={record['energy_mev']:g} MeV "
                f"angle={record['angle_deg']:g} deg runtime={record['runtime_seconds']:.2f} s",
                flush=True,
            )

    records.sort(key=lambda row: (row["particle"], row["energy_mev"], row["angle_deg"]))
    manifest = {
        "protocol": "PENTAGON-GEANT4-PARTICLE-RESPONSE-GRID-V1",
        "driver": str(driver),
        "driver_sha256": file_sha256(driver),
        "python": sys.version,
        "events_em": args.events_em,
        "events_other": args.events_other,
        "angles_deg": list(ANGLES_DEG),
        "em_energies_mev": list(EM_ENERGIES_MEV),
        "other_energies_mev": list(OTHER_ENERGIES_MEV),
        "started_unix": started,
        "finished_unix": time.time(),
        "records": records,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "GEANT4_PARTICLE_RESPONSE_GRID_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    failures = [record for record in records if record["status"] == "failed"]
    print(f"manifest={manifest_path} failures={len(failures)}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
