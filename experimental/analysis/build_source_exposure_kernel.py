#!/usr/bin/env python3
"""Build the actual-live-time zenith exposure kernel for catalog sources.

The kernel is independent of the detector-response simulation.  It records how
many live seconds each fixed ICRS direction spent in fine true-zenith bins over
the 44 V12 DAQ runs.  A later calculation can therefore fold any simulated
effective-area curve through the real observing history without approximating
the archive as a uniform calendar year.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd


# These are the exact constants used in the definitive V12 reconstruction.
SITE_LATITUDE_DEG = 35.0 + 43.0 / 60.0
SITE_LONGITUDE_DEG = 51.0 + 20.0 / 60.0
TEHRAN_ZONE = ZoneInfo("Asia/Tehran")
J2000_UNIX_NS = 946_728_000_000_000_000.0
DAY_NS = 86_400_000_000_000.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def days_since_j2000_from_ns(utc_ns: np.ndarray) -> np.ndarray:
    return (np.asarray(utc_ns, dtype=np.float64) - J2000_UNIX_NS) / DAY_NS


def gmst_degrees_from_ns(utc_ns: np.ndarray) -> np.ndarray:
    days = days_since_j2000_from_ns(utc_ns)
    centuries = days / 36525.0
    return (
        280.46061837
        + 360.98564736629 * days
        + 0.000387933 * centuries**2
        - centuries**3 / 38_710_000.0
    ) % 360.0


def precession_angles_from_ns(utc_ns: np.ndarray) -> tuple[np.ndarray, ...]:
    centuries = days_since_j2000_from_ns(utc_ns) / 36525.0
    scale = math.pi / (180.0 * 3600.0)
    zeta = (
        2306.2181 * centuries
        + 0.30188 * centuries**2
        + 0.017998 * centuries**3
    ) * scale
    z_angle = (
        2306.2181 * centuries
        + 1.09468 * centuries**2
        + 0.018203 * centuries**3
    ) * scale
    theta = (
        2004.3109 * centuries
        - 0.42665 * centuries**2
        - 0.041833 * centuries**3
    ) * scale
    return zeta, z_angle, theta


def icrs_to_mean_of_date_vectors(
    ra_deg: np.ndarray, dec_deg: np.ndarray, utc_ns: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Invert the exact mean-of-date-to-ICRS rotation implemented in V12."""

    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))[None, :]
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))[None, :]
    x0 = np.cos(dec) * np.cos(ra)
    y0 = np.cos(dec) * np.sin(ra)
    z1 = np.sin(dec)

    zeta, z_angle, theta = precession_angles_from_ns(utc_ns)
    czeta, szeta = np.cos(zeta)[:, None], np.sin(zeta)[:, None]
    ctheta, stheta = np.cos(theta)[:, None], np.sin(theta)[:, None]
    cz, sz = np.cos(z_angle)[:, None], np.sin(z_angle)[:, None]

    # Inverse of Rz(-zeta) Ry(theta) Rz(-z).
    x1 = czeta * x0 - szeta * y0
    y2 = szeta * x0 + czeta * y0
    x2 = ctheta * x1 - stheta * z1
    z = stheta * x1 + ctheta * z1
    x = cz * x2 - sz * y2
    y = sz * x2 + cz * y2
    return np.arctan2(y, x), np.arcsin(np.clip(z, -1.0, 1.0))


def midpoint_samples(
    start_local: str, duration_seconds: float, step_seconds: float
) -> tuple[np.ndarray, np.ndarray]:
    start = pd.Timestamp(start_local).tz_localize(
        TEHRAN_ZONE, ambiguous="raise", nonexistent="raise"
    )
    start_ns = int(start.tz_convert("UTC").value)
    edges = np.arange(0.0, duration_seconds, step_seconds, dtype=float)
    widths = np.minimum(step_seconds, duration_seconds - edges)
    midpoints_ns = start_ns + np.rint((edges + 0.5 * widths) * 1e9).astype(np.int64)
    return midpoints_ns, widths


def zenith_angles_deg(
    ra_deg: np.ndarray, dec_deg: np.ndarray, utc_ns: np.ndarray
) -> np.ndarray:
    ra_mod, dec_mod = icrs_to_mean_of_date_vectors(ra_deg, dec_deg, utc_ns)
    lst = np.deg2rad(
        (gmst_degrees_from_ns(utc_ns) + SITE_LONGITUDE_DEG) % 360.0
    )[:, None]
    hour_angle = lst - ra_mod
    latitude = math.radians(SITE_LATITUDE_DEG)
    cosine = (
        math.sin(latitude) * np.sin(dec_mod)
        + math.cos(latitude) * np.cos(dec_mod) * np.cos(hour_angle)
    )
    return np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-summary",
        type=Path,
        default=Path(__file__).parents[1] / "audit" / "v12" / "run_summary.csv",
    )
    parser.add_argument(
        "--sources",
        type=Path,
        default=Path(__file__).parents[1] / "audit" / "lhaaso_37_sources.csv",
    )
    parser.add_argument("--time-step-seconds", type=float, default=30.0)
    parser.add_argument("--zenith-bin-deg", type=float, default=0.25)
    parser.add_argument("--maximum-zenith-deg", type=float, default=45.0)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).with_name("source_exposure"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.time_step_seconds <= 0 or args.zenith_bin_deg <= 0:
        raise ValueError("time and zenith steps must be positive")
    if args.maximum_zenith_deg <= 0 or args.maximum_zenith_deg > 90:
        raise ValueError("maximum zenith must lie in (0, 90]")

    run_path = args.run_summary.resolve()
    source_path = args.sources.resolve()
    runs = pd.read_csv(run_path)
    sources = pd.read_csv(source_path)
    required_run = {"run", "start_local", "duration_hours"}
    required_source = {"source", "ra_deg", "dec_deg"}
    if not required_run.issubset(runs.columns):
        raise ValueError(f"run summary missing {sorted(required_run - set(runs.columns))}")
    if not required_source.issubset(sources.columns):
        raise ValueError(f"source table missing {sorted(required_source - set(sources.columns))}")

    edges = np.arange(
        0.0,
        args.maximum_zenith_deg + 0.5 * args.zenith_bin_deg,
        args.zenith_bin_deg,
    )
    if not np.isclose(edges[-1], args.maximum_zenith_deg):
        raise ValueError("maximum zenith must be divisible by zenith bin width")
    hist_seconds = np.zeros((len(sources), len(edges) - 1), dtype=float)
    all_live_seconds = 0.0
    sampled_seconds = 0.0
    sample_count = 0

    ra = sources.ra_deg.to_numpy(float)
    dec = sources.dec_deg.to_numpy(float)
    for run in runs.itertuples(index=False):
        duration_seconds = float(run.duration_hours) * 3600.0
        times_ns, widths = midpoint_samples(
            str(run.start_local), duration_seconds, args.time_step_seconds
        )
        theta = zenith_angles_deg(ra, dec, times_ns)
        for source_index in range(len(sources)):
            hist_seconds[source_index] += np.histogram(
                theta[:, source_index], bins=edges, weights=widths
            )[0]
        all_live_seconds += duration_seconds
        sampled_seconds += float(widths.sum())
        sample_count += len(widths)

    if not np.isclose(all_live_seconds, sampled_seconds, rtol=0, atol=1e-6):
        raise RuntimeError("midpoint quadrature failed to preserve DAQ live time")

    centers = 0.5 * (edges[:-1] + edges[1:])
    long_records = []
    summary_records = []
    for source_index, source in sources.iterrows():
        exposure = hist_seconds[source_index]
        visible = float(exposure.sum())
        for bin_index, seconds in enumerate(exposure):
            long_records.append(
                {
                    "source": source.source,
                    "ra_deg": float(source.ra_deg),
                    "dec_deg": float(source.dec_deg),
                    "zenith_low_deg": float(edges[bin_index]),
                    "zenith_high_deg": float(edges[bin_index + 1]),
                    "zenith_center_deg": float(centers[bin_index]),
                    "live_seconds": float(seconds),
                }
            )
        summary_records.append(
            {
                "source": source.source,
                "ra_deg": float(source.ra_deg),
                "dec_deg": float(source.dec_deg),
                "archive_live_hours": all_live_seconds / 3600.0,
                "visible_live_hours_theta_lt_45": visible / 3600.0,
                "visible_fraction_of_archive_live_time": visible / all_live_seconds,
                "exposure_weighted_mean_zenith_deg": float(
                    np.dot(exposure, centers) / visible
                )
                if visible
                else np.nan,
            }
        )

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    kernel_path = output_dir / "lhaaso_37_actual_live_zenith_kernel.csv"
    summary_path = output_dir / "lhaaso_37_actual_live_exposure_summary.csv"
    manifest_path = output_dir / "SOURCE_EXPOSURE_KERNEL_MANIFEST.json"
    pd.DataFrame(long_records).to_csv(kernel_path, index=False)
    summary_table = pd.DataFrame(summary_records).sort_values(
        "visible_live_hours_theta_lt_45", ascending=False
    )
    summary_table.to_csv(summary_path, index=False)
    manifest = {
        "protocol": "PENTAGON-SOURCE-EXPOSURE-KERNEL-V1",
        "run_summary": str(run_path),
        "run_summary_sha256": sha256(run_path),
        "sources": str(source_path),
        "sources_sha256": sha256(source_path),
        "site_latitude_deg": SITE_LATITUDE_DEG,
        "site_longitude_deg": SITE_LONGITUDE_DEG,
        "timezone": str(TEHRAN_ZONE),
        "coordinate_method": (
            "inverse of the V12 mean-of-date-to-ICRS precession rotation, "
            "followed by the V12 GMST expression"
        ),
        "time_quadrature": "midpoint",
        "time_step_seconds": args.time_step_seconds,
        "zenith_bin_deg": args.zenith_bin_deg,
        "maximum_zenith_deg": args.maximum_zenith_deg,
        "run_count": int(len(runs)),
        "source_count": int(len(sources)),
        "sample_count": sample_count,
        "archive_live_seconds": all_live_seconds,
        "kernel_csv": str(kernel_path),
        "kernel_csv_sha256": sha256(kernel_path),
        "summary_csv": str(summary_path),
        "summary_csv_sha256": sha256(summary_path),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(summary_table.to_string(index=False))
    print(f"manifest={manifest_path}")


if __name__ == "__main__":
    main()
