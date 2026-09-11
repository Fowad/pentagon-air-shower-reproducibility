#!/usr/bin/env python3
"""Exact V12 tests of externally reported 2017 Cygnus X-3 gamma-ray windows."""

from pathlib import Path

import pandas as pd

from exact_permutation_moments import load_v12, ROOT


WINDOWS = (
    (
        "Fermi-LAT MJD 57798-57801",
        "2017-02-14T00:00:00Z",
        3.0,
        "ATel 10109; nearly 7 sigma above 100 MeV",
    ),
    (
        "AGILE 2017-02-27 to 2017-03-01",
        "2017-02-27T03:00:00Z",
        2.0,
        "ATel 10138; near 4 sigma above 100 MeV",
    ),
    (
        "AGILE 2017-03-15",
        "2017-03-15T00:00:00Z",
        1.0,
        "ATel 10179; near 4 sigma above 100 MeV",
    ),
    (
        "Fermi-LAT 2017-04-03",
        "2017-04-03T00:00:00Z",
        1.0,
        "ATel 10243; gamma-ray activity at onset of major radio flare",
    ),
)


def main():
    v12 = load_v12()
    events = pd.read_csv(ROOT / "audit" / "v12" / "reconstructed_events.csv.gz")
    config = v12.Config(
        data_dir=Path("."),
        output_dir=Path("."),
        event_cache=None,
        run_summary_cache=None,
        targeted_source_realizations=10_000,
        make_figures=False,
    )
    rows = []
    for offset, (label, start, duration, context) in enumerate(WINDOWS):
        print(label, flush=True)
        result = v12._exact_window_target_statistic(
            events,
            config,
            308.1166,
            40.9434,
            start,
            duration,
            15_000 + offset,
        )
        rows.append(
            {
                "window": label,
                "start_utc": start,
                "duration_days": duration,
                "external_context": context,
                **result,
            }
        )
    out = pd.DataFrame(rows)
    out["bonferroni_four_windows"] = (4 * out["one_sided_p"]).clip(upper=1.0)
    destination = ROOT / "audit" / "results" / "cygx3_external_windows_exact.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(destination, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
