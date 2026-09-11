#!/usr/bin/env python3
"""Fold one unthinned CORSIKA shower through the five-detector response model."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation, shift, uniform_filter
from scipy.spatial import cKDTree
import yaml


C_M_PER_S = 299_792_458.0
MIP_MPV_MEV = 3.57852
THRESHOLDS_MIP = np.asarray([0.25, 0.5, 1.0, 2.0])
SOURCE_APERTURE_RADIUS_DEG = 8.0
DETECTOR_COORDINATES_M = np.asarray(
    [
        [-2.7586, -3.6756],
        [-4.3133, 1.0717],
        [-0.2780, 4.0196],
        [3.7707, 1.0942],
        [2.2376, -3.6617],
    ],
    dtype=float,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def particle_class(pdg: int) -> str | None:
    if pdg == 22:
        return "gamma"
    if pdg == 11:
        return "e-"
    if pdg == -11:
        return "e+"
    if abs(pdg) == 13:
        return "mu-"
    if abs(pdg) == 2212:
        return "proton"
    if abs(pdg) == 2112:
        return "neutron"
    if abs(pdg) in (211, 321):
        return "pi+"
    return None


class ResponseLibrary:
    def __init__(
        self,
        grid_dir: Path,
        bare_gamma_grid: Path | None = None,
        low_energy_extension: Path | None = None,
        low_energy_other_extension: Path | None = None,
        high_energy_em_extension: Path | None = None,
    ):
        def response_path(record: dict, directory: Path) -> Path:
            """Resolve a manifest output both in-place and after archive relocation."""
            recorded = Path(record["output"])
            if recorded.exists():
                return recorded
            relocated = directory / recorded.name
            if relocated.exists():
                return relocated
            raise FileNotFoundError(
                f"GEANT4 response file not found at {recorded} or {relocated}"
            )

        self.grid_dir = grid_dir.resolve()
        manifest_path = self.grid_dir / "GEANT4_PARTICLE_RESPONSE_GRID_MANIFEST.json"
        self.manifest_sha256 = sha256(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        if len(manifest["records"]) != 270:
            raise RuntimeError("GEANT4 response grid is incomplete")
        self.samples: dict[tuple[str, float, float], np.ndarray] = {}
        for record in manifest["records"]:
            path = response_path(record, self.grid_dir)
            if sha256(path) != record["sha256"]:
                raise RuntimeError(f"GEANT4 grid hash mismatch: {path}")
            with np.load(path, allow_pickle=False) as data:
                key = (
                    str(data["particle"]),
                    float(data["kinetic_energy_mev"]),
                    float(data["angle_deg"]),
                )
                if not bool(data["include_steel"]):
                    raise RuntimeError(f"expected steel-enclosed response: {path}")
                self.samples[key] = np.asarray(data["deposited_energy_mev"], dtype=float)
        self.gamma_geometry = "steel_enclosure"
        self.bare_gamma_manifest_sha256 = None
        self.low_energy_extension_manifest_sha256 = None
        self.low_energy_other_extension_manifest_sha256 = None
        self.high_energy_em_extension_manifest_sha256 = None
        if bare_gamma_grid is not None:
            bare_gamma_grid = bare_gamma_grid.resolve()
            bare_manifest_path = bare_gamma_grid / "GEANT4_BARE_GAMMA_GRID_MANIFEST.json"
            self.bare_gamma_manifest_sha256 = sha256(bare_manifest_path)
            bare_manifest = json.loads(bare_manifest_path.read_text())
            if len(bare_manifest["records"]) != 50:
                raise RuntimeError("bare-gamma response grid is incomplete")
            for record in bare_manifest["records"]:
                path = response_path(record, bare_gamma_grid)
                if sha256(path) != record["sha256"]:
                    raise RuntimeError(f"bare-gamma grid hash mismatch: {path}")
                with np.load(path, allow_pickle=False) as data:
                    if str(data["particle"]) != "gamma" or bool(data["include_steel"]):
                        raise RuntimeError(f"invalid bare-gamma response: {path}")
                    key = (
                        "gamma",
                        float(data["kinetic_energy_mev"]),
                        float(data["angle_deg"]),
                    )
                    self.samples[key] = np.asarray(
                        data["deposited_energy_mev"], dtype=float
                    )
            self.gamma_geometry = "bare_plastic"
        if low_energy_other_extension is not None:
            low_energy_other_extension = low_energy_other_extension.resolve()
            other_manifest_path = (
                low_energy_other_extension
                / "GEANT4_LOW_ENERGY_OTHER_EXTENSION_MANIFEST.json"
            )
            self.low_energy_other_extension_manifest_sha256 = sha256(
                other_manifest_path
            )
            extension = json.loads(other_manifest_path.read_text())
            if len(extension["records"]) != 80:
                raise RuntimeError("low-energy muon/hadron extension is incomplete")
            for record in extension["records"]:
                path = response_path(record, low_energy_other_extension)
                if sha256(path) != record["sha256"]:
                    raise RuntimeError(f"low-energy other-grid hash mismatch: {path}")
                with np.load(path, allow_pickle=False) as data:
                    key = (
                        str(data["particle"]),
                        float(data["kinetic_energy_mev"]),
                        float(data["angle_deg"]),
                    )
                    if not bool(data["include_steel"]):
                        raise RuntimeError(f"low-energy other grid lacks steel: {path}")
                    self.samples[key] = np.asarray(
                        data["deposited_energy_mev"], dtype=float
                    )
        if low_energy_extension is not None:
            low_energy_extension = low_energy_extension.resolve()
            extension_manifest_path = (
                low_energy_extension / "GEANT4_LOW_ENERGY_EM_EXTENSION_MANIFEST.json"
            )
            self.low_energy_extension_manifest_sha256 = sha256(
                extension_manifest_path
            )
            extension = json.loads(extension_manifest_path.read_text())
            if len(extension["records"]) != 60:
                raise RuntimeError("low-energy GEANT4 extension is incomplete")
            for record in extension["records"]:
                particle = str(record["particle"])
                bare = bool(record["bare"])
                wanted = (particle != "gamma" and not bare) or (
                    particle == "gamma"
                    and bare == (self.gamma_geometry == "bare_plastic")
                )
                if not wanted:
                    continue
                path = response_path(record, low_energy_extension)
                if sha256(path) != record["sha256"]:
                    raise RuntimeError(f"low-energy grid hash mismatch: {path}")
                with np.load(path, allow_pickle=False) as data:
                    key = (
                        str(data["particle"]),
                        float(data["kinetic_energy_mev"]),
                        float(data["angle_deg"]),
                    )
                    if bool(data["include_steel"]) == bare:
                        raise RuntimeError(f"low-energy geometry mismatch: {path}")
                    self.samples[key] = np.asarray(
                        data["deposited_energy_mev"], dtype=float
                    )
        if high_energy_em_extension is not None:
            high_energy_em_extension = high_energy_em_extension.resolve()
            high_manifest_path = (
                high_energy_em_extension
                / "GEANT4_HIGH_ENERGY_EM_EXTENSION_MANIFEST.json"
            )
            self.high_energy_em_extension_manifest_sha256 = sha256(
                high_manifest_path
            )
            extension = json.loads(high_manifest_path.read_text())
            if len(extension["records"]) != 40:
                raise RuntimeError("high-energy electromagnetic extension is incomplete")
            for record in extension["records"]:
                particle = str(record["particle"])
                bare = bool(record["bare"])
                wanted = (particle != "gamma" and not bare) or (
                    particle == "gamma"
                    and bare == (self.gamma_geometry == "bare_plastic")
                )
                if not wanted:
                    continue
                path = response_path(record, high_energy_em_extension)
                if sha256(path) != record["sha256"]:
                    raise RuntimeError(f"high-energy grid hash mismatch: {path}")
                with np.load(path, allow_pickle=False) as data:
                    key = (
                        str(data["particle"]),
                        float(data["kinetic_energy_mev"]),
                        float(data["angle_deg"]),
                    )
                    if bool(data["include_steel"]) == bare:
                        raise RuntimeError(f"high-energy geometry mismatch: {path}")
                    self.samples[key] = np.asarray(
                        data["deposited_energy_mev"], dtype=float
                    )
        self.axes: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for particle in sorted({key[0] for key in self.samples}):
            energies = np.asarray(sorted({key[1] for key in self.samples if key[0] == particle}))
            angles = np.asarray(sorted({key[2] for key in self.samples if key[0] == particle}))
            self.axes[particle] = energies, angles

    @staticmethod
    def _mixture_indices(
        values: np.ndarray, grid: np.ndarray, transform, rng: np.random.Generator
    ) -> tuple[np.ndarray, int, int]:
        transformed = transform(values)
        grid_transformed = transform(grid)
        upper = np.searchsorted(grid_transformed, transformed, side="right")
        below = int(np.count_nonzero(upper == 0))
        above = int(np.count_nonzero(upper == len(grid)))
        upper = np.clip(upper, 1, len(grid) - 1)
        lower = upper - 1
        denominator = grid_transformed[upper] - grid_transformed[lower]
        fraction = np.divide(
            transformed - grid_transformed[lower],
            denominator,
            out=np.zeros_like(transformed, dtype=float),
            where=denominator != 0,
        )
        choose_upper = rng.random(len(values)) < np.clip(fraction, 0, 1)
        indices = np.where(choose_upper, upper, lower)
        indices[values <= grid[0]] = 0
        indices[values >= grid[-1]] = len(grid) - 1
        return indices, below, above

    def sample(
        self,
        pdg: np.ndarray,
        kinetic_energy_gev: np.ndarray,
        nz: np.ndarray,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, dict]:
        deposits = np.zeros(len(pdg), dtype=float)
        classes = np.asarray([particle_class(int(value)) for value in pdg], dtype=object)
        audit = {
            "unsupported_particles": int(np.count_nonzero(classes == None)),  # noqa: E711
            "upward_particles": int(np.count_nonzero(nz >= 0)),
            "by_class": {},
        }
        downward_angle = np.rad2deg(np.arccos(np.clip(-nz, 0, 1)))
        energy_mev = kinetic_energy_gev * 1_000.0
        for particle, (energy_grid, angle_grid) in self.axes.items():
            selected = np.flatnonzero(classes == particle)
            if not len(selected):
                continue
            e_index, e_below, e_above = self._mixture_indices(
                energy_mev[selected], energy_grid, np.log, rng
            )
            secant = lambda value: 1.0 / np.cos(np.deg2rad(np.clip(value, 0, 89.9)))
            a_index, a_below, a_above = self._mixture_indices(
                downward_angle[selected], angle_grid, secant, rng
            )
            for ei in range(len(energy_grid)):
                for ai in range(len(angle_grid)):
                    local = np.flatnonzero((e_index == ei) & (a_index == ai))
                    if not len(local):
                        continue
                    target = selected[local]
                    samples = self.samples[(particle, float(energy_grid[ei]), float(angle_grid[ai]))]
                    deposits[target] = samples[rng.integers(0, len(samples), size=len(target))]
            audit["by_class"][particle] = {
                "particles": int(len(selected)),
                "energy_below_grid": e_below,
                "energy_above_grid": e_above,
                "angle_below_grid": a_below,
                "angle_above_grid": a_above,
                "nonzero_deposits": int(np.count_nonzero(deposits[selected] > 0)),
            }
        deposits[nz >= 0] = 0.0
        return deposits, audit


def prompt_trigger_times(
    times: np.ndarray,
    deposits: np.ndarray,
    thresholds_mev: np.ndarray,
    window_seconds: float,
) -> np.ndarray:
    output = np.full(len(thresholds_mev), np.nan, dtype=float)
    if not len(times):
        return output
    order = np.argsort(times)
    t = times[order]
    e = deposits[order]
    left = 0
    running = 0.0
    for right in range(len(t)):
        running += e[right]
        while t[right] - t[left] > window_seconds:
            running -= e[left]
            left += 1
        newly_crossed = np.isnan(output) & (running >= thresholds_mev)
        output[newly_crossed] = t[right]
        if np.all(np.isfinite(output)):
            break
    return output


def reconstruct_directions(trigger_times: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # trigger_times: candidate x station
    baselines = DETECTOR_COORDINATES_M[1:] - DETECTOR_COORDINATES_M[0]
    inverse = np.linalg.pinv(baselines)
    relative_distance = (
        trigger_times[:, 1:] - trigger_times[:, [0]]
    ) * C_M_PER_S
    horizontal = relative_distance @ inverse.T
    norm = np.linalg.norm(horizontal, axis=1)
    valid = np.isfinite(norm) & (norm <= 1)
    theta = np.full(len(norm), np.nan)
    theta[valid] = np.rad2deg(np.arcsin(norm[valid]))
    return horizontal, theta, valid


def fast_candidate_centers(
    x: np.ndarray,
    y: np.ndarray,
    deposits: np.ndarray,
    extent_m: float,
    step_m: float,
) -> tuple[np.ndarray, dict]:
    bins = int(round(2 * extent_m / step_m))
    if not np.isclose(bins * step_m, 2 * extent_m):
        raise ValueError("2*extent must be divisible by the core-grid spacing")
    detector_bins = int(round(0.5 / step_m))
    if detector_bins < 1 or not np.isclose(detector_bins * step_m, 0.5):
        raise ValueError("0.5 m detector width must be divisible by grid spacing")
    histogram, _, _ = np.histogram2d(
        y,
        x,
        bins=bins,
        range=((-extent_m, extent_m), (-extent_m, extent_m)),
        weights=deposits,
    )
    detector_deposit = (
        uniform_filter(histogram, size=detector_bins, mode="constant")
        * detector_bins
        * detector_bins
    )
    # A deliberately loose prefilter is followed by exact square-footprint and
    # pulse-time calculations.  Dilation prevents a bin edge from excluding a
    # physically valid core point.
    prefilter_mev = 0.5 * THRESHOLDS_MIP[0] * MIP_MPV_MEV
    candidate = np.ones_like(detector_deposit, dtype=bool)
    for dx, dy in DETECTOR_COORDINATES_M:
        station_map = shift(
            detector_deposit,
            shift=(-dy / step_m, -dx / step_m),
            order=1,
            mode="constant",
            cval=0,
            prefilter=False,
        )
        station_possible = binary_dilation(station_map >= prefilter_mev, iterations=2)
        candidate &= station_possible
    iy, ix = np.nonzero(candidate)
    centers = np.column_stack(
        [
            -extent_m + (ix + 0.5) * step_m,
            -extent_m + (iy + 0.5) * step_m,
        ]
    )
    audit = {
        "histogram_bins": bins,
        "detector_bins": detector_bins,
        "prefilter_mev": prefilter_mev,
        "candidate_core_points": int(len(centers)),
        "nonzero_particle_deposits": int(np.count_nonzero(deposits > 0)),
        "deposited_energy_in_grid_mev": float(histogram.sum()),
        "candidate_touches_grid_boundary": bool(
            len(ix)
            and (
                np.min(ix) <= 2
                or np.max(ix) >= bins - 3
                or np.min(iy) <= 2
                or np.max(iy) >= bins - 3
            )
        ),
    }
    return centers, audit


def evaluate_replica(
    particles: pd.DataFrame,
    response: ResponseLibrary,
    response_seed: int,
    extent_m: float,
    step_m: float,
    pulse_windows_ns: tuple[float, ...],
    electronics_sigma_ns: float,
    primary_direction: np.ndarray,
    query_workers: int,
) -> tuple[list[dict], dict]:
    rng = np.random.default_rng(response_seed)
    deposits, sampling_audit = response.sample(
        particles.pdg.to_numpy(int),
        particles.kinetic_energy.to_numpy(float),
        particles.nz.to_numpy(float),
        rng,
    )
    x = particles.x.to_numpy(float)
    y = particles.y.to_numpy(float)
    times = particles.time.to_numpy(float)
    centers, candidate_audit = fast_candidate_centers(
        x, y, deposits, extent_m, step_m
    )
    active = deposits > 0
    tree = cKDTree(np.column_stack([x[active], y[active]]))
    active_times = times[active]
    active_deposits = deposits[active]
    thresholds_mev = THRESHOLDS_MIP * MIP_MPV_MEV
    windows = np.asarray(pulse_windows_ns) * 1e-9
    trigger_times = np.full(
        (len(windows), len(centers), len(DETECTOR_COORDINATES_M), len(thresholds_mev)),
        np.nan,
        dtype=float,
    )
    for station, detector_coordinate in enumerate(DETECTOR_COORDINATES_M):
        detector_centers = centers + detector_coordinate
        neighborhoods = tree.query_ball_point(
            detector_centers, r=0.25, p=np.inf, workers=query_workers
        )
        for core_index, indices in enumerate(neighborhoods):
            if not len(indices):
                continue
            index = np.asarray(indices, dtype=int)
            for window_index, window in enumerate(windows):
                trigger_times[window_index, core_index, station] = prompt_trigger_times(
                    active_times[index], active_deposits[index], thresholds_mev, window
                )

    records: list[dict] = []
    truth_horizontal = primary_direction[:2]
    truth_cosine = -primary_direction[2]
    for window_index, window_ns in enumerate(pulse_windows_ns):
        for threshold_index, threshold_mip in enumerate(THRESHOLDS_MIP):
            station_times = trigger_times[window_index, :, :, threshold_index]
            fivefold = np.all(np.isfinite(station_times), axis=1)
            coincidence_span = np.full(len(centers), np.nan)
            coincidence_span[fivefold] = np.ptp(station_times[fivefold], axis=1)
            coincidence = fivefold & (coincidence_span <= 200e-9)
            selected_indices = np.flatnonzero(coincidence)
            theta = np.full(len(centers), np.nan)
            separation = np.full(len(centers), np.nan)
            direction_valid = np.zeros(len(centers), dtype=bool)
            if len(selected_indices):
                jittered = station_times[selected_indices] + rng.normal(
                    0, electronics_sigma_ns * 1e-9, size=(len(selected_indices), 5)
                )
                horizontal, reconstructed_theta, valid = reconstruct_directions(jittered)
                theta[selected_indices] = reconstructed_theta
                direction_valid[selected_indices] = valid
                valid_indices = selected_indices[valid]
                if len(valid_indices):
                    reconstructed_horizontal = horizontal[valid]
                    reconstructed_vertical = np.sqrt(
                        np.maximum(0, 1 - np.sum(reconstructed_horizontal**2, axis=1))
                    )
                    cosine = (
                        reconstructed_horizontal @ truth_horizontal
                        + reconstructed_vertical * truth_cosine
                    )
                    separation[valid_indices] = np.rad2deg(
                        np.arccos(np.clip(cosine, -1, 1))
                    )
            accepted = coincidence & direction_valid & (theta < 45)
            source_aperture = accepted & (separation <= SOURCE_APERTURE_RADIUS_DEG)
            accepted_separation = separation[accepted]
            records.append(
                {
                    "response_seed": response_seed,
                    "pulse_window_ns": float(window_ns),
                    "threshold_mip": float(threshold_mip),
                    "candidate_core_points": int(len(centers)),
                    "fivefold_core_points": int(np.count_nonzero(fivefold)),
                    "coincidence_core_points": int(np.count_nonzero(coincidence)),
                    "accepted_core_points": int(np.count_nonzero(accepted)),
                    "effective_area_m2": float(np.count_nonzero(accepted) * step_m**2),
                    "source_aperture_core_points": int(
                        np.count_nonzero(source_aperture)
                    ),
                    "source_aperture_effective_area_m2": float(
                        np.count_nonzero(source_aperture) * step_m**2
                    ),
                    "source_aperture_containment_fraction": float(
                        np.count_nonzero(source_aperture) / np.count_nonzero(accepted)
                    )
                    if np.count_nonzero(accepted)
                    else np.nan,
                    "coincidence_rejection_fraction": float(
                        1 - np.count_nonzero(coincidence) / np.count_nonzero(fivefold)
                    )
                    if np.count_nonzero(fivefold)
                    else np.nan,
                    "direction_rejection_fraction": float(
                        1 - np.count_nonzero(accepted) / np.count_nonzero(coincidence)
                    )
                    if np.count_nonzero(coincidence)
                    else np.nan,
                    "median_angular_error_deg": float(np.nanmedian(accepted_separation))
                    if len(accepted_separation)
                    else np.nan,
                    "p68_angular_error_deg": float(np.nanquantile(accepted_separation, 0.68))
                    if len(accepted_separation)
                    else np.nan,
                    "max_accepted_core_radius_m": float(
                        np.max(np.linalg.norm(centers[accepted], axis=1))
                    )
                    if np.any(accepted)
                    else np.nan,
                }
            )
    audit = {
        "sampling": sampling_audit,
        "candidate_grid": candidate_audit,
        "deposited_energy_total_mev": float(deposits.sum()),
        "deposit_nonzero_fraction": float(np.mean(deposits > 0)),
    }
    return records, audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shower-dir", type=Path, required=True)
    parser.add_argument(
        "--response-grid",
        type=Path,
        default=Path(__file__).with_name("geant4_response") / "particle_grid_v1",
    )
    parser.add_argument(
        "--bare-gamma-grid",
        type=Path,
        help="replace only the photon response with the bare-plastic systematic grid",
    )
    parser.add_argument(
        "--low-energy-extension",
        type=Path,
        help="optional GEANT4 0.05--0.30 MeV electromagnetic response extension",
    )
    parser.add_argument(
        "--low-energy-other-extension",
        type=Path,
        help="optional GEANT4 20--200 MeV muon/hadron response extension",
    )
    parser.add_argument(
        "--high-energy-em-extension",
        type=Path,
        help="optional GEANT4 30--100 GeV electromagnetic response extension",
    )
    parser.add_argument("--extent-m", type=float, default=150.0)
    parser.add_argument("--step-m", type=float, default=0.25)
    parser.add_argument("--replicas", type=int, default=3)
    parser.add_argument("--base-seed", type=int, default=2026082001)
    parser.add_argument("--pulse-windows-ns", type=float, nargs="+", default=[5, 10, 20])
    parser.add_argument("--electronics-sigma-ns", type=float, default=1.85)
    parser.add_argument(
        "--query-workers",
        type=int,
        default=1,
        help="workers used by each KD-tree neighborhood query",
    )
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.query_workers <= 0:
        raise ValueError("query-workers must be positive")
    shower_dir = args.shower_dir.resolve()
    particles_path = shower_dir / "particles" / "particles.parquet"
    particles = pd.read_parquet(particles_path)
    if len(particles) == 0:
        raise RuntimeError("CORSIKA shower has no ground particles")
    if not np.allclose(particles.weight.to_numpy(float), 1.0):
        raise RuntimeError("this program accepts only unthinned, unit-weight showers")
    primary = yaml.safe_load((shower_dir / "primary" / "summary.yaml").read_text())[
        "shower_0"
    ]
    primary_direction = np.asarray(
        [primary["nx"], primary["ny"], primary["nz"]], dtype=float
    )
    response = ResponseLibrary(
        args.response_grid,
        args.bare_gamma_grid,
        args.low_energy_extension,
        args.low_energy_other_extension,
        args.high_energy_em_extension,
    )
    records: list[dict] = []
    audits: list[dict] = []
    for replica in range(args.replicas):
        seed = args.base_seed + replica
        replica_records, audit = evaluate_replica(
            particles,
            response,
            seed,
            args.extent_m,
            args.step_m,
            tuple(args.pulse_windows_ns),
            args.electronics_sigma_ns,
            primary_direction,
            args.query_workers,
        )
        records.extend(replica_records)
        audits.append({"response_seed": seed, **audit})
        print(f"completed response replica {replica + 1}/{args.replicas}", flush=True)

    output_prefix = args.output_prefix.resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_prefix.with_suffix(".csv")
    json_path = output_prefix.with_suffix(".json")
    table = pd.DataFrame(records).sort_values(
        ["response_seed", "pulse_window_ns", "threshold_mip"]
    )
    table.to_csv(csv_path, index=False)
    summary = {
        "protocol": "PENTAGON-UNTHINNED-EFFECTIVE-AREA-V1",
        "shower_dir": str(shower_dir),
        "particles_path": str(particles_path),
        "particles_sha256": sha256(particles_path),
        "ground_particles": int(len(particles)),
        "all_particle_weights_one": True,
        "primary": primary,
        "response_grid": str(args.response_grid.resolve()),
        "response_grid_manifest_sha256": response.manifest_sha256,
        "gamma_geometry": response.gamma_geometry,
        "bare_gamma_manifest_sha256": response.bare_gamma_manifest_sha256,
        "low_energy_extension_manifest_sha256": (
            response.low_energy_extension_manifest_sha256
        ),
        "low_energy_other_extension_manifest_sha256": (
            response.low_energy_other_extension_manifest_sha256
        ),
        "high_energy_em_extension_manifest_sha256": (
            response.high_energy_em_extension_manifest_sha256
        ),
        "mip_mpv_mev": MIP_MPV_MEV,
        "source_aperture_radius_deg": SOURCE_APERTURE_RADIUS_DEG,
        "thresholds_mip": THRESHOLDS_MIP.tolist(),
        "pulse_windows_ns": list(args.pulse_windows_ns),
        "electronics_sigma_ns": args.electronics_sigma_ns,
        "query_workers": args.query_workers,
        "extent_m": args.extent_m,
        "step_m": args.step_m,
        "replicas": args.replicas,
        "records_csv": str(csv_path),
        "records_csv_sha256": sha256(csv_path),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "audits": audits,
    }
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(table.to_string(index=False))
    print(f"summary={json_path}")


if __name__ == "__main__":
    main()
