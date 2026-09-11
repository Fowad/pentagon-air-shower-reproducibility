#!/usr/bin/env python3
"""GEANT4 deposited-energy response of one pentagon-array scintillator.

The model is deliberately detector-local.  It propagates one specified ground
particle through a 50 x 50 x 2 cm PVT scintillator enclosed by 1 mm iron.  The
result is the distribution of deposited energy in the active plastic,
conditional on the particle trajectory crossing the detector's horizontal
mid-plane inside its nominal 50 x 50 cm footprint.

This program does not model the unknown discriminator, PMT collection, or
five-fold trigger.  Its output is intended for the explicit threshold and
efficiency envelope defined in UHE_ANALYSIS_PREREGISTRATION.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

from geant4_pybind import (
    FTFP_BERT,
    G4Box,
    G4LogicalVolume,
    G4NistManager,
    G4PVPlacement,
    G4ParticleGun,
    G4ParticleTable,
    G4Random,
    G4RunManager,
    G4RunManagerFactory,
    G4RunManagerType,
    G4ThreeVector,
    G4UImanager,
    G4UserEventAction,
    G4UserSteppingAction,
    G4VERSION_NUMBER,
    G4VUserDetectorConstruction,
    G4VUserPrimaryGeneratorAction,
    MTwistEngine,
    MeV,
    cm,
    deg,
    mm,
)


PARTICLE_ALIASES = {
    "gamma": "gamma",
    "photon": "gamma",
    "e-": "e-",
    "electron": "e-",
    "e+": "e+",
    "positron": "e+",
    "mu-": "mu-",
    "mu+": "mu+",
    "proton": "proton",
    "neutron": "neutron",
    "pi+": "pi+",
    "pi-": "pi-",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class PentagonDetectorConstruction(G4VUserDetectorConstruction):
    """Active PVT slab plus the documented 1 mm galvanized enclosure."""

    def __init__(self, include_steel: bool = True):
        super().__init__()
        self.include_steel = include_steel
        self.scoring_volume = None

    def Construct(self):
        nist = G4NistManager.Instance()
        air = nist.FindOrBuildMaterial("G4_AIR")
        pvt = nist.FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE")
        iron = nist.FindOrBuildMaterial("G4_Fe")

        world_solid = G4Box("World", 1.5 * 100 * cm, 1.5 * 100 * cm, 1.5 * 100 * cm)
        world_logic = G4LogicalVolume(world_solid, air, "World")
        world_physical = G4PVPlacement(
            None, G4ThreeVector(), world_logic, "World", None, False, 0, False
        )

        scintillator_solid = G4Box("Scintillator", 25 * cm, 25 * cm, 1 * cm)
        scintillator_logic = G4LogicalVolume(
            scintillator_solid, pvt, "Scintillator"
        )
        G4PVPlacement(
            None,
            G4ThreeVector(),
            scintillator_logic,
            "Scintillator",
            world_logic,
            False,
            0,
            False,
        )
        self.scoring_volume = scintillator_logic

        if self.include_steel:
            half_skin = 0.5 * mm
            plate_xy = G4Box("SteelTopBottom", 25.1 * cm, 25.1 * cm, half_skin)
            plate_logic = G4LogicalVolume(plate_xy, iron, "SteelTopBottom")
            for copy_number, z in enumerate((1 * cm + half_skin, -1 * cm - half_skin)):
                G4PVPlacement(
                    None,
                    G4ThreeVector(0, 0, z),
                    plate_logic,
                    "SteelTopBottom",
                    world_logic,
                    False,
                    copy_number,
                    False,
                )

            side_x = G4Box("SteelSideX", half_skin, 25 * cm, 1 * cm)
            side_x_logic = G4LogicalVolume(side_x, iron, "SteelSideX")
            for copy_number, x in enumerate((25 * cm + half_skin, -25 * cm - half_skin)):
                G4PVPlacement(
                    None,
                    G4ThreeVector(x, 0, 0),
                    side_x_logic,
                    "SteelSideX",
                    world_logic,
                    False,
                    copy_number,
                    False,
                )

            side_y = G4Box("SteelSideY", 25 * cm, half_skin, 1 * cm)
            side_y_logic = G4LogicalVolume(side_y, iron, "SteelSideY")
            for copy_number, y in enumerate((25 * cm + half_skin, -25 * cm - half_skin)):
                G4PVPlacement(
                    None,
                    G4ThreeVector(0, y, 0),
                    side_y_logic,
                    "SteelSideY",
                    world_logic,
                    False,
                    copy_number,
                    False,
                )

        return world_physical


class ParticleGenerator(G4VUserPrimaryGeneratorAction):
    def __init__(self, particle: str, energy_mev: float, angle_deg: float):
        super().__init__()
        self.gun = G4ParticleGun(1)
        definition = G4ParticleTable.GetParticleTable().FindParticle(particle)
        if definition is None:
            raise ValueError(f"GEANT4 particle is unavailable: {particle}")
        self.gun.SetParticleDefinition(definition)
        self.gun.SetParticleEnergy(energy_mev * MeV)

        angle = angle_deg * deg
        self.nx = float(np.sin(angle_deg * np.pi / 180.0))
        self.nz = -float(np.cos(angle_deg * np.pi / 180.0))
        self.gun.SetParticleMomentumDirection(G4ThreeVector(self.nx, 0, self.nz))
        self.start_z = 30 * cm

    def GeneratePrimaries(self, event):
        # Uniform mid-plane crossing point.  Moving the launch point upstream
        # preserves that crossing coordinate for inclined trajectories.
        x_cross = (2.0 * np.random.random() - 1.0) * 25 * cm
        y_cross = (2.0 * np.random.random() - 1.0) * 25 * cm
        # With r(s) = r_start + n*s, the mid-plane is reached at
        # s = -start_z/nz.  Therefore x_start must be displaced upstream by
        # +(nx/nz)*start_z so that x(s) is exactly x_cross at z=0.
        x_start = x_cross + (self.nx / self.nz) * self.start_z
        self.gun.SetParticlePosition(G4ThreeVector(x_start, y_cross, self.start_z))
        self.gun.GeneratePrimaryVertex(event)


class DepositEventAction(G4UserEventAction):
    def __init__(self):
        super().__init__()
        self.current = 0.0
        self.deposits_mev: list[float] = []

    def BeginOfEventAction(self, event):
        self.current = 0.0

    def EndOfEventAction(self, event):
        self.deposits_mev.append(float(self.current / MeV))


class DepositSteppingAction(G4UserSteppingAction):
    def __init__(self, detector: PentagonDetectorConstruction, event_action: DepositEventAction):
        super().__init__()
        self.detector = detector
        self.event_action = event_action

    def UserSteppingAction(self, step):
        pre_volume = step.GetPreStepPoint().GetTouchable().GetVolume()
        if pre_volume is None:
            return
        if pre_volume.GetLogicalVolume() == self.detector.scoring_volume:
            self.event_action.current += step.GetTotalEnergyDeposit()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particle", default="mu-", choices=sorted(PARTICLE_ALIASES))
    parser.add_argument("--energy-mev", type=float, default=10_000.0)
    parser.add_argument("--angle-deg", type=float, default=0.0)
    parser.add_argument("--events", type=int, default=1_000)
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--bare", action="store_true", help="omit the 1 mm steel enclosure")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.energy_mev <= 0 or args.events <= 0:
        raise ValueError("energy and event count must be positive")
    if not 0 <= args.angle_deg < 80:
        raise ValueError("angle must be in [0, 80) degrees")

    particle = PARTICLE_ALIASES[args.particle]
    np.random.seed(args.seed)
    # Keep the Python owner alive for the full run; GEANT4 stores a pointer to
    # the engine rather than making a copy.
    random_engine = MTwistEngine()
    G4Random.setTheEngine(random_engine)
    G4Random.setTheSeed(args.seed)

    run_manager = G4RunManagerFactory.CreateRunManager(G4RunManagerType.Serial)
    detector = PentagonDetectorConstruction(include_steel=not args.bare)
    physics = FTFP_BERT()
    physics.SetVerboseLevel(0)
    physics.SetDefaultCutValue(0.1 * mm)
    run_manager.SetUserInitialization(detector)
    run_manager.SetUserInitialization(physics)

    # GEANT4 requires the physics list to be assigned before any user-action
    # object is instantiated.
    generator = ParticleGenerator(particle, args.energy_mev, args.angle_deg)
    event_action = DepositEventAction()
    stepping_action = DepositSteppingAction(detector, event_action)

    run_manager.SetUserAction(generator)
    run_manager.SetUserAction(event_action)
    run_manager.SetUserAction(stepping_action)
    run_manager.Initialize()

    ui = G4UImanager.GetUIpointer()
    for command in ("/run/verbose 0", "/event/verbose 0", "/tracking/verbose 0"):
        ui.ApplyCommand(command)
    run_manager.BeamOn(args.events)

    deposits = np.asarray(event_action.deposits_mev, dtype=np.float64)
    if len(deposits) != args.events:
        raise RuntimeError(f"GEANT4 returned {len(deposits)} deposits for {args.events} events")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        deposited_energy_mev=deposits,
        particle=np.asarray(particle),
        kinetic_energy_mev=np.asarray(args.energy_mev),
        angle_deg=np.asarray(args.angle_deg),
        seed=np.asarray(args.seed),
        include_steel=np.asarray(not args.bare),
        geant4_version_number=np.asarray(G4VERSION_NUMBER),
    )
    quantiles = np.quantile(deposits, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1])
    summary = {
        "particle": particle,
        "kinetic_energy_mev": args.energy_mev,
        "angle_deg": args.angle_deg,
        "events": args.events,
        "seed": args.seed,
        "include_steel": not args.bare,
        "geant4_version_number": int(G4VERSION_NUMBER),
        "nonzero_fraction": float(np.mean(deposits > 0)),
        "mean_deposit_mev": float(np.mean(deposits)),
        "standard_deviation_mev": float(np.std(deposits, ddof=1)),
        "quantile_probabilities": [0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1],
        "quantiles_mev": quantiles.tolist(),
        "output_sha256": sha256(args.output),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    sys.stdout.flush()
    sys.stderr.flush()
    # geant4_pybind 0.1.3 can fault while Python and GEANT4 independently
    # dismantle user-action objects at interpreter shutdown.  All run output
    # is finalized above; bypassing that binding-only teardown is deterministic
    # and does not skip any event or file finalization.
    os._exit(0)


if __name__ == "__main__":
    main()
