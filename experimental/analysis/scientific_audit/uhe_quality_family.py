#!/usr/bin/env python3
"""Locked reconstruction-quality UHE source-family test.

Implements Test 1 in UHE_ANALYSIS_PREREGISTRATION.md.  The calculation uses
exact first and second moments of the V12 blockwise permutation distribution;
it does not use Monte Carlo estimates for an individual test's background.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
V12_PATH = (
    ROOT
    / "audit"
    / "handoff"
    / "PENTAGON_V12_NEW_WORK_HANDOFF_PACKAGE"
    / "Pentagon_Array_Directional_Analysis_V12.py"
)
EVENTS_PATH = ROOT / "audit" / "v12" / "reconstructed_events.csv.gz"
RUNS_PATH = ROOT / "audit" / "v12" / "run_summary.csv"
SOURCES_PATH = ROOT / "audit" / "pevatron_sources_working.csv"
PRIMARY_PATH = ROOT / "audit" / "v12" / "results" / "primary_results_V12.json"
OUTPUT_DIR = ROOT / "investigation" / "results"

EXPECTED_SHA256 = {
    EVENTS_PATH: "2ca83067b8b14ee0c1a82eb8b4c51eef63c4b1d2b6e0526c9ddaa7ef6f8ebb40",
    RUNS_PATH: "384c46f9b9a30f3a853b22054ec1337acc06e876bce96b9401da5acafa361a25",
    SOURCES_PATH: "a43e824ad832bbd375148d6f186b7d798ef4e78e00c51dcb381a3eeb93fc0fee",
}
QUALITY_LABELS = ("all", "best_50pct", "best_25pct", "best_10pct")
QUALITY_FRACTIONS = (1.0, 0.50, 0.25, 0.10)
PERIOD_LABELS = ("full",) + tuple(f"interval_{i}" for i in range(1, 7))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inputs() -> None:
    for path, expected in EXPECTED_SHA256.items():
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"Input identity mismatch for {path}: {observed}")


def load_v12():
    spec = importlib.util.spec_from_file_location("pentagon_v12_uhe_quality", V12_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {V12_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def source_vectors(sources: pd.DataFrame) -> np.ndarray:
    ra = np.radians(sources["ra_deg"].to_numpy(float))
    dec = np.radians(sources["dec_deg"].to_numpy(float))
    return np.column_stack(
        (np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))
    )


def plane_residual_rms_ns(events: pd.DataFrame, runs: pd.DataFrame) -> np.ndarray:
    centers = runs.set_index("run")
    center_matrix = np.column_stack(
        [
            events["run"].map(centers[f"center_{channel}"]).to_numpy(float)
            for channel in ("ch21", "ch31", "ch41", "ch51")
        ]
    )
    channels = events[["ch21", "ch31", "ch41", "ch51"]].to_numpy(float)
    delays_ns = (channels - center_matrix) * (200.0 / 1024.0)
    detector_coordinates = np.array(
        [
            [-2.7586, -3.6756],
            [-4.3133, 1.0717],
            [-0.2780, 4.0196],
            [3.7707, 1.0942],
            [2.2376, -3.6617],
        ],
        dtype=float,
    )
    geometry = detector_coordinates[1:] - detector_coordinates[0]
    c_m_per_ns = 299_702_547.0e-9
    horizontal_solution = c_m_per_ns * delays_ns @ np.linalg.pinv(geometry).T
    predicted_ns = horizontal_solution @ geometry.T / c_m_per_ns
    residual_ns = delays_ns - predicted_ns
    return np.sqrt(np.mean(np.square(residual_ns), axis=1))


def runwise_quality_masks(
    events: pd.DataFrame, residual_rms_ns: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame]:
    masks = np.zeros((len(events), len(QUALITY_FRACTIONS)), dtype=bool)
    masks[:, 0] = True
    rows = []
    for run, index in events.groupby("run", sort=False).indices.items():
        index = np.asarray(index, dtype=np.int64)
        order = index[np.argsort(residual_rms_ns[index], kind="stable")]
        record = {
            "run": run,
            "accepted_events": int(len(index)),
            "residual_median_ns": float(np.median(residual_rms_ns[index])),
        }
        for column, (label, fraction) in enumerate(
            zip(QUALITY_LABELS[1:], QUALITY_FRACTIONS[1:]), start=1
        ):
            count = int(math.floor(fraction * len(index)))
            masks[order[:count], column] = True
            record[f"{label}_events"] = count
            record[f"{label}_edge_ns"] = float(residual_rms_ns[order[count - 1]])
        rows.append(record)
    if not np.all(masks[:, 1] >= masks[:, 2]) or not np.all(masks[:, 2] >= masks[:, 3]):
        raise RuntimeError("Quality masks are not nested")
    return masks, pd.DataFrame(rows)


def observed_counts(
    v12,
    events: pd.DataFrame,
    sources: pd.DataFrame,
    interval_id: np.ndarray,
    quality_masks: np.ndarray,
) -> np.ndarray:
    # Shape: quality, period, source.
    out = np.zeros((len(QUALITY_LABELS), len(PERIOD_LABELS), len(sources)), dtype=np.int64)
    ra = events["ra_deg"].to_numpy(float)
    dec = events["dec_deg"].to_numpy(float)
    for source_index, source in enumerate(sources.itertuples(index=False)):
        aperture = v12._exact_circle_membership(
            ra, dec, source.ra_deg, source.dec_deg, 8.0
        )
        for quality_index in range(len(QUALITY_LABELS)):
            selected = aperture & quality_masks[:, quality_index]
            out[quality_index, 0, source_index] = np.count_nonzero(selected)
            for interval in range(6):
                out[quality_index, interval + 1, source_index] = np.count_nonzero(
                    selected & (interval_id == interval)
                )
    return out


def exact_moments(v12, state, sources: pd.DataFrame, quality_masks: np.ndarray):
    # Shape: quality, period, source.
    shape = (len(QUALITY_LABELS), len(PERIOD_LABELS), len(sources))
    means = np.zeros(shape, dtype=np.float64)
    variances = np.zeros(shape, dtype=np.float64)
    vectors = source_vectors(sources)
    cos_radius = math.cos(math.radians(8.0)) - 1.0e-15

    for block_index, indices in enumerate(state.blocks):
        indices = np.asarray(indices, dtype=np.int64)
        n = len(indices)
        hour_angle = state.hour_angle_deg[indices][:, None]
        dec_mean = np.broadcast_to(state.declination_mean_deg[indices][:, None], (n, n))
        lst = np.broadcast_to(state.observed_lst_deg[indices][None, :], (n, n))
        angles = [
            np.broadcast_to(angle[indices][None, :], (n, n))
            for angle in state.precession
        ]
        ra, dec = v12.mean_of_date_to_icrs_with_angles(
            (lst - hour_angle) % 360.0, dec_mean, *angles
        )
        ra_rad = np.radians(ra.ravel())
        dec_rad = np.radians(dec.ravel())
        sky_vectors = np.column_stack(
            (
                np.cos(dec_rad) * np.cos(ra_rad),
                np.cos(dec_rad) * np.sin(ra_rad),
                np.sin(dec_rad),
            )
        )
        membership = ((sky_vectors @ vectors.T) >= cos_radius).reshape(
            n, n, len(sources)
        )
        period_masks = [np.ones(n, dtype=bool)] + [
            state.interval_id[indices] == interval for interval in range(6)
        ]

        for quality_index in range(len(QUALITY_LABELS)):
            row_mask = quality_masks[indices, quality_index]
            if not np.any(row_mask):
                continue
            for period_index, column_mask in enumerate(period_masks):
                if not np.any(column_mask):
                    continue
                selected = membership[:, column_mask, :] & row_mask[:, None, None]
                total = selected.sum(axis=(0, 1), dtype=np.int64).astype(float)
                row_sum = selected.sum(axis=1, dtype=np.int64).astype(float)
                column_sum = selected.sum(axis=0, dtype=np.int64).astype(float)
                means[quality_index, period_index] += total / n
                if n > 1:
                    centered_ss = (
                        total
                        - np.square(row_sum).sum(axis=0) / n
                        - np.square(column_sum).sum(axis=0) / n
                        + np.square(total) / (n * n)
                    )
                    variances[quality_index, period_index] += centered_ss / (n - 1)

        if block_index % 100 == 0 or block_index + 1 == len(state.blocks):
            print(f"blocks {block_index + 1}/{len(state.blocks)}", flush=True)

    return means, variances


def main() -> None:
    verify_inputs()
    v12 = load_v12()
    events = pd.read_csv(EVENTS_PATH)
    runs = pd.read_csv(RUNS_PATH)
    sources = pd.read_csv(SOURCES_PATH)
    sources = sources[sources["p300_core"].eq("yes")].reset_index(drop=True)
    if len(sources) != 17:
        raise RuntimeError(f"Locked source family has {len(sources)} rows, expected 17")

    primary = json.loads(PRIMARY_PATH.read_text())
    boundaries = pd.to_datetime(primary["boundaries"], utc=True, format="mixed")
    times = pd.to_datetime(events["utc_time"], utc=True, format="mixed")
    interval_id = np.searchsorted(
        boundaries.asi8[1:-1], times.astype("int64"), side="right"
    ).astype(np.int16)

    residual = plane_residual_rms_ns(events, runs)
    quality_masks, quality_audit = runwise_quality_masks(events, residual)
    config = v12.Config(
        data_dir=Path("."),
        output_dir=Path("."),
        event_cache=None,
        run_summary_cache=None,
        make_figures=False,
    )
    state = v12.randomization_state(events, config, interval_id)
    observed = observed_counts(v12, events, sources, interval_id, quality_masks)
    means, variances = exact_moments(v12, state, sources, quality_masks)
    standard_deviations = np.sqrt(np.maximum(variances, 0.0))
    z = np.divide(
        observed - means,
        standard_deviations,
        out=np.full_like(means, np.nan),
        where=standard_deviations > 0,
    )
    local_p = norm.sf(z)

    rows = []
    for quality_index, quality_label in enumerate(QUALITY_LABELS):
        for period_index, period_label in enumerate(PERIOD_LABELS):
            for source_index, source in sources.iterrows():
                rows.append(
                    {
                        "quality_sample": quality_label,
                        "period": period_label,
                        "source": source["source"],
                        "ra_deg": source["ra_deg"],
                        "dec_deg": source["dec_deg"],
                        "N": int(observed[quality_index, period_index, source_index]),
                        "B_exact": means[quality_index, period_index, source_index],
                        "sB_exact": standard_deviations[
                            quality_index, period_index, source_index
                        ],
                        "Z_exact": z[quality_index, period_index, source_index],
                        "normal_one_sided_p": local_p[
                            quality_index, period_index, source_index
                        ],
                    }
                )
    table = pd.DataFrame(rows).sort_values("Z_exact", ascending=False)
    number_of_tests = int(np.count_nonzero(np.isfinite(z)))
    maximum = table.iloc[0]
    minimum_local_p = float(maximum["normal_one_sided_p"])
    bonferroni = min(1.0, number_of_tests * minimum_local_p)
    followup_required = bool(float(maximum["Z_exact"]) >= 3.5 or bonferroni <= 0.25)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT_DIR / "uhe_quality_source_family_statistics.csv", index=False)
    quality_audit.to_csv(OUTPUT_DIR / "plane_residual_quality_audit_by_run.csv", index=False)
    np.savez_compressed(
        OUTPUT_DIR / "uhe_quality_source_family_moments.npz",
        observed=observed,
        mean=means,
        variance=variances,
        z=z,
        local_p=local_p,
        residual_rms_ns=residual,
        quality_masks=quality_masks,
        quality_labels=np.asarray(QUALITY_LABELS),
        period_labels=np.asarray(PERIOD_LABELS),
        source_names=sources["source"].to_numpy(str),
    )
    summary = {
        "locked_test_count": 17 * 7 * 4,
        "finite_test_count": number_of_tests,
        "maximum": {
            key: (value.item() if hasattr(value, "item") else value)
            for key, value in maximum.to_dict().items()
        },
        "bonferroni_family_upper_bound": bonferroni,
        "correlation_aware_followup_required_by_locked_rule": followup_required,
        "residual_rms_ns_quantiles": {
            str(q): float(np.quantile(residual, q))
            for q in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
        },
        "selected_event_counts": {
            label: int(quality_masks[:, index].sum())
            for index, label in enumerate(QUALITY_LABELS)
        },
        "input_sha256": {str(path.relative_to(ROOT)): digest for path, digest in EXPECTED_SHA256.items()},
    }
    (OUTPUT_DIR / "uhe_quality_source_family_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    print("\nTop ten tests:\n")
    print(table.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
