from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta


ROOT = Path(__file__).resolve().parent
V12 = ROOT / "v12"
RESULTS = V12 / "results"
HANDOFF = ROOT / "handoff" / "PENTAGON_V12_NEW_WORK_HANDOFF_PACKAGE"


def plus_one_p(null: np.ndarray, threshold: float) -> tuple[int, int, float]:
    values = np.asarray(null, dtype=float)
    finite = np.isfinite(values)
    k = int(finite.sum())
    exceedances = int(np.sum(values[finite] >= float(threshold)))
    return exceedances, k, float((exceedances + 1) / (k + 1))


def binomial_interval(exceedances: int, trials: int, alpha: float = 0.05) -> tuple[float, float]:
    # Exact Clopper-Pearson interval for the underlying exceedance probability.
    lo = 0.0 if exceedances == 0 else float(beta.ppf(alpha / 2, exceedances, trials - exceedances + 1))
    hi = 1.0 if exceedances == trials else float(beta.ppf(1 - alpha / 2, exceedances + 1, trials - exceedances))
    return lo, hi


def event_fingerprint(events: pd.DataFrame, run_order: list[str]) -> str:
    digest = hashlib.sha256()
    columns = ["daq_counter", "ch21", "ch31", "ch41", "ch51"]
    run_values = events["run"].astype(str)
    for name in run_order:
        subset = events.loc[run_values.eq(name)]
        digest.update(name.encode("ascii") + b"\0")
        digest.update(np.ascontiguousarray(subset[columns].to_numpy(dtype="<i8")).tobytes())
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads((V12 / "V12_RUN_MANIFEST.json").read_text())
    primary = json.loads((RESULTS / "primary_results_V12.json").read_text())
    peaks = pd.read_csv(RESULTS / "directional_excess_peaks_V12.csv")
    run_summary = pd.read_csv(V12 / "run_summary.csv")

    code = HANDOFF / "Pentagon_Array_Directional_Analysis_V12.py"
    code_bytes = code.stat().st_size
    code_hash = hashlib.sha256(code.read_bytes()).hexdigest()

    events = pd.read_csv(
        V12 / "reconstructed_events.csv.gz",
        usecols=["run", "daq_counter", "ch21", "ch31", "ch41", "ch51"],
        dtype={
            "run": "string",
            "daq_counter": "int64",
            "ch21": "int64",
            "ch31": "int64",
            "ch41": "int64",
            "ch51": "int64",
        },
    )
    fingerprint = event_fingerprint(events, run_summary["run"].astype(str).tolist())

    with np.load(V12 / "mc_checkpoint_reference_V12.npz", allow_pickle=False) as ref:
        ref_meta = {key: ref[key].item() for key in [
            "state_fingerprint", "operator_fingerprint", "code_sha256", "release_id",
            "aperture_radius_deg", "target_realizations", "completed_realizations",
            "random_seed", "nside",
        ]}
        means = ref["means"].copy()
        std = np.sqrt(ref["m2"].copy() / (int(ref_meta["completed_realizations"]) - 1))
        cumulative_means = ref["cumulative_means"].copy()
        cumulative_std = np.sqrt(
            ref["cumulative_m2"].copy() / (int(ref_meta["completed_realizations"]) - 1)
        )

    with np.load(V12 / "mc_checkpoint_calibration_V12.npz", allow_pickle=False) as cal:
        cal_meta = {key: cal[key].item() for key in [
            "state_fingerprint", "analysis_fingerprint", "code_sha256", "release_id",
            "target_realizations", "completed_realizations", "random_seed", "nside",
        ]}
        arrays = {key: cal[key].copy() for key in [
            "null_full_max", "null_local_primary", "null_interval_max",
            "null_cumulative_primary", "null_gamma_full_all", "null_gamma_interval_all",
            "null_gamma_full_high", "null_gamma_interval_high",
        ]}

    checkpoint_array_health = {
        name: {
            "shape": list(value.shape),
            "finite": int(np.isfinite(value).sum()),
            "total": int(value.size),
            "min": float(np.nanmin(value)),
            "max": float(np.nanmax(value)),
            "mean": float(np.nanmean(value)),
            "std": float(np.nanstd(value, ddof=1)),
        }
        for name, value in arrays.items()
    }

    peak_crosschecks = []
    for row_index, row in peaks.reset_index(drop=True).iterrows():
        pixel = int(row.pixel)
        recomputed_b = float(means[row_index, pixel])
        recomputed_s = float(std[row_index, pixel])
        recomputed_z = float((int(row.n_p) - recomputed_b) / recomputed_s)
        peak_crosschecks.append({
            "period": row.period,
            "pixel": pixel,
            "delta_B": recomputed_b - float(row.b_p),
            "delta_sB": recomputed_s - float(row.reference_std),
            "delta_Z": recomputed_z - float(row.z),
        })

    observed_full = float(peaks.loc[peaks.period.eq("full"), "z"].iloc[0])
    observed_intervals = peaks.loc[peaks.period.ne("full"), "z"].to_numpy(float)
    p_recomputed: dict[str, object] = {}
    for name, values, threshold in [
        ("primary_local_p", arrays["null_local_primary"], observed_full),
        ("primary_map_p", arrays["null_full_max"], observed_full),
        ("six_interval_family_p", np.max(arrays["null_interval_max"], axis=1), float(observed_intervals.max())),
    ]:
        exceed, k, p = plus_one_p(values, threshold)
        p_recomputed[name] = {
            "threshold": threshold,
            "exceedances": exceed,
            "trials": k,
            "p": p,
            "reported": float(primary[name]),
            "delta": p - float(primary[name]),
            "mc_95pct_interval": list(binomial_interval(exceed, k)),
        }

    interval_map = []
    for i, threshold in enumerate(observed_intervals):
        exceed, k, p = plus_one_p(arrays["null_interval_max"][:, i], threshold)
        interval_map.append({
            "interval": i + 1,
            "threshold": float(threshold),
            "exceedances": exceed,
            "p": p,
            "reported": float(peaks.loc[peaks.period.eq(f"interval_{i+1}"), "map_p"].iloc[0]),
            "mc_95pct_interval": list(binomial_interval(exceed, k)),
        })
    p_recomputed["interval_map_p"] = interval_map

    family_specs = [
        (
            "gamma_catalog_all_valid",
            arrays["null_gamma_full_all"],
            arrays["null_gamma_interval_all"],
            np.asarray(primary["gamma_catalog_all_valid_observed_max"], float),
        ),
        (
            "gamma_catalog_high_occupancy",
            arrays["null_gamma_full_high"],
            arrays["null_gamma_interval_high"],
            np.asarray(primary["gamma_catalog_high_occupancy_observed_max"], float),
        ),
    ]
    family_recomputed = {}
    for label, nf, ni, obs in family_specs:
        full_exceed, k, full_p = plus_one_p(nf, obs[0])
        interval_values = []
        for i in range(6):
            exceed, _, p = plus_one_p(ni[:, i], obs[i + 1])
            interval_values.append({
                "interval": i + 1,
                "threshold": float(obs[i + 1]),
                "exceedances": exceed,
                "p": p,
                "mc_95pct_interval": list(binomial_interval(exceed, k)),
            })
        six_exceed, _, six_p = plus_one_p(np.nanmax(ni, axis=1), float(np.nanmax(obs[1:])))
        family_recomputed[label] = {
            "full": {
                "threshold": float(obs[0]),
                "exceedances": full_exceed,
                "p": full_p,
                "mc_95pct_interval": list(binomial_interval(full_exceed, k)),
            },
            "intervals": interval_values,
            "six_interval": {
                "threshold": float(np.nanmax(obs[1:])),
                "exceedances": six_exceed,
                "p": six_p,
                "mc_95pct_interval": list(binomial_interval(six_exceed, k)),
            },
        }

    old_threshold_checks = {
        "full_old_z_on_v12_null": plus_one_p(arrays["null_full_max"], 2.598345),
        "interval3_old_z_on_v12_null": plus_one_p(arrays["null_interval_max"][:, 2], 3.968750),
        "six_interval_old_z_on_v12_null": plus_one_p(np.max(arrays["null_interval_max"], axis=1), 3.968750),
    }

    report = {
        "identity": {
            "manifest_stage": manifest["stage"],
            "code_bytes_actual": code_bytes,
            "code_bytes_manifest": manifest["code_bytes"],
            "code_sha256_actual": code_hash,
            "code_sha256_manifest": manifest["code_sha256"],
            "release_manifest": manifest["release"],
            "accepted_events_actual": int(len(events)),
            "accepted_events_manifest": int(manifest["accepted_events"]),
            "run_count_actual": int(run_summary.shape[0]),
            "run_count_unique_events": int(events["run"].nunique()),
            "livetime_sum_hours": float(run_summary["duration_hours"].sum()),
            "livetime_manifest_hours": float(manifest["livetime_hours"]),
            "event_sha256_actual": fingerprint,
            "event_sha256_manifest": manifest["event_sha256"],
            "accepted_counts_sum": int(run_summary["accepted_event_count"].sum()),
        },
        "checkpoint_metadata": {"reference": ref_meta, "calibration": cal_meta},
        "checkpoint_array_health": checkpoint_array_health,
        "reference_array_health": {
            "means_shape": list(means.shape),
            "std_shape": list(std.shape),
            "cumulative_means_shape": list(cumulative_means.shape),
            "cumulative_std_shape": list(cumulative_std.shape),
            "means_all_finite": bool(np.all(np.isfinite(means))),
            "std_all_finite": bool(np.all(np.isfinite(std))),
            "std_min": float(std.min()),
            "std_max": float(std.max()),
        },
        "peak_crosschecks": peak_crosschecks,
        "p_recomputed": p_recomputed,
        "catalog_families_recomputed": family_recomputed,
        "old_thresholds_evaluated_on_v12_null": old_threshold_checks,
    }
    print(json.dumps(report, indent=2, default=lambda x: x.item() if hasattr(x, "item") else x))


if __name__ == "__main__":
    main()
