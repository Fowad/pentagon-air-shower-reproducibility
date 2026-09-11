#!/usr/bin/env python3
"""Direct 10k V12 permutation check of the LS I +61 303 phase-window lead."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from exact_permutation_moments import ROOT, load_v12
from lsi_orbital_phase_exact import (
    ORBIT_DAYS,
    PHASE_EPOCH_MJD,
    SOURCE_DEC_DEG,
    SOURCE_RA_DEG,
)


def main():
    v12 = load_v12()
    events = pd.read_csv(ROOT / "audit" / "v12" / "reconstructed_events.csv.gz")
    times = pd.to_datetime(events["utc_time"], utc=True, format="mixed")
    mjd = times.astype("int64").to_numpy(float) / 1.0e9 / 86400.0 + 40587.0
    phase = np.mod(mjd - PHASE_EPOCH_MJD, ORBIT_DAYS) / ORBIT_DAYS
    target = (phase >= 0.3) & (phase < 0.6)
    config = v12.Config(
        data_dir=Path("."),
        output_dir=Path("."),
        event_cache=None,
        run_summary_cache=None,
        targeted_source_realizations=10_000,
        make_figures=False,
    )
    state = v12.randomization_state(events, config, np.zeros(len(events), dtype=np.int16))
    relevant = []
    for block_index, indices in enumerate(state.blocks):
        column_target = target[indices]
        if not np.any(column_target):
            continue
        matrix = v12._exact_target_block_matrix(
            state,
            indices,
            SOURCE_RA_DEG,
            SOURCE_DEC_DEG,
            config.aperture_radius_deg,
        )
        relevant.append((matrix, column_target, np.arange(len(indices), dtype=np.int64)))
        if block_index % 100 == 0:
            print(f"prepared block {block_index + 1}/{len(state.blocks)}", flush=True)
    observed = int(
        np.count_nonzero(
            v12._exact_circle_membership(
                events.loc[target, "ra_deg"].to_numpy(float),
                events.loc[target, "dec_deg"].to_numpy(float),
                SOURCE_RA_DEG,
                SOURCE_DEC_DEG,
                config.aperture_radius_deg,
            )
        )
    )
    values = np.empty(config.targeted_source_realizations, dtype=np.int32)
    for trial in range(config.targeted_source_realizations):
        rng = v12.realization_rng(config.random_seed, 61303, trial)
        count = 0
        for matrix, column_target, row_index in relevant:
            permutation = rng.permutation(len(row_index))
            count += int(
                np.sum(matrix[row_index, permutation] & column_target[permutation])
            )
        values[trial] = count
        if (trial + 1) % 1000 == 0:
            print(f"trials {trial + 1}/{config.targeted_source_realizations}", flush=True)
    background = float(values.mean())
    background_std = float(values.std(ddof=1))
    result = {
        "source": "LS I +61 303",
        "test": "published 25-100 TeV phase 0.3-0.6",
        "trials": config.targeted_source_realizations,
        "N": observed,
        "B": background,
        "sB": background_std,
        "Z": float((observed - background) / background_std),
        "one_sided_p": float(
            (1 + np.count_nonzero(values >= observed))
            / (config.targeted_source_realizations + 1)
        ),
    }
    destination = ROOT / "audit" / "results" / "lsi_phase_0p3_0p6_direct_mc.json"
    destination.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
