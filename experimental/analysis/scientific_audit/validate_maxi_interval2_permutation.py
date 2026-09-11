#!/usr/bin/env python3
"""Direct 10k V12 permutation validation of the leading P300 local statistic."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from exact_permutation_moments import ROOT, load_v12


def main():
    v12 = load_v12()
    events = pd.read_csv(ROOT / "audit" / "v12" / "reconstructed_events.csv.gz")
    primary = json.loads(
        (ROOT / "audit" / "v12" / "results" / "primary_results_V12.json").read_text()
    )
    boundaries = pd.to_datetime(primary["boundaries"], utc=True, format="mixed")
    times = pd.to_datetime(events["utc_time"], utc=True, format="mixed")
    interval_id = np.searchsorted(
        boundaries.asi8[1:-1], times.astype("int64"), side="right"
    ).astype(np.int16)
    config = v12.Config(
        data_dir=Path("."),
        output_dir=Path("."),
        event_cache=None,
        run_summary_cache=None,
        targeted_source_realizations=10_000,
        make_figures=False,
    )
    state = v12.randomization_state(events, config, interval_id)
    n, b, sb, z, p = v12._exact_interval_target_statistic(
        events,
        state,
        interval_id,
        1,
        275.22,
        7.39,
        config,
        16_000,
    )
    result = {
        "source": "MAXI J1820+070 UHE component",
        "period": "interval_2",
        "trials": 10_000,
        "N": n,
        "B": b,
        "sB": sb,
        "Z": z,
        "one_sided_p": p,
    }
    destination = ROOT / "audit" / "results" / "maxi_j1820_interval2_direct_mc.json"
    destination.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
