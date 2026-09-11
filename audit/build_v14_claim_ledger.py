#!/usr/bin/env python3
"""Promote the V13 claim ledger to V14 using the newly supplied evidence."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "prior_v13_evidence" / "V13_MASTER_NUMERICAL_CLAIM_LEDGER.csv"
OUTPUT = ROOT / "fresh_evidence" / "V14_MASTER_NUMERICAL_CLAIM_LEDGER.csv"


RESOLVED = {
    "ABS-006": (
        "Independent 1,000,000-draw correlated-Gaussian rerun p=0.717741; frozen 5,000,000-draw value p=0.717866",
        "Resolved by source_family_probability_recheck.json; both values round to 0.718.",
    ),
    "ABS-007": (
        "17 directions and 102 interval cells; independent p=0.138106; frozen p=0.138331",
        "Resolved from exact covariance with a new random stream.",
    ),
    "REF-001": (
        "The 2 h, approximately 4 h, and 8 h examples remain supported; V14 adds the ARGO-YBJ citation",
        "The V13 citation omission is corrected in V14.",
    ),
    "SIM-001": (
        "All 18 R3 logs record CORSIKA8 0.1.99 and the declared models/site; fresh detector run records GEANT4 version 1141",
        "Resolved by compact-log audit and fresh GEANT4 run.",
    ),
    "SIM-004": (
        "Fresh generator records [50,50,2] cm plastic dimensions and 1 mm iron enclosure",
        "Resolved by fresh 790,000-event GEANT4 basis and event metadata.",
    ),
    "SIM-005": (
        "Fresh 20,000-muon MIP fit gives 3.577887617 MeV versus archived 3.57852 MeV",
        "Difference is -0.000632383 MeV (-0.0177%); package tolerance passes.",
    ),
    "LHA-002": (
        "Independent full-family p=0.717741; equal stack p=0.673640; 100-TeV-flux stack p=0.705442",
        "Resolved from exact covariance and independently regenerated Gaussian family trials.",
    ),
    "UHE-001": (
        "17 directions; max Z=1.914951; independent full-family p=0.355444",
        "Resolved from exact covariance with a new 1,000,000-draw stream.",
    ),
    "UHE-002": (
        "102 cells; max Z=2.962586; independent family p=0.138106",
        "Resolved from exact covariance with a new 1,000,000-draw stream.",
    ),
    "UHE-003": (
        "Independent equal-source six-interval p=0.130726; aggregate-count p=0.151250",
        "Both reproduce the frozen 0.131--0.151 range.",
    ),
    "SRC-001": (
        "LS I interval-1 Z=2.962586; direct p=0.00157998; controlling 17x6 p=0.138106",
        "Local excess is real as a fluctuation but not significant after its parent family.",
    ),
    "SRC-002": (
        "MAXI J1820+070 interval-2 Z=2.496898; direct p=0.00714993; parent-family p=0.138331",
        "Resolved by the definitive 100,000-trial checkpoint and summary.",
    ),
    "ORB-001": (
        "LS I direct local p=0.00157998 (standalone), 0.00179998 (joint grid); source-six p=0.00974990",
        "Resolved by direct 100,000-trial arrays with add-one counting.",
    ),
    "ORB-002": (
        "Ten phase bins p=0.0930591; two windows p=0.0319497; 12 cells p=0.0738393; 60 cells p=0.115389",
        "Resolved by direct 100,000-trial arrays with add-one counting.",
    ),
    "CON-003": (
        "All named source-family and phase-family probabilities independently reproduce and remain non-significant at their controlling family level",
        "The V14 null conclusion is supported.",
    ),
}


def main() -> None:
    frame = pd.read_csv(SOURCE)
    # These two statements were deleted in V14 rather than carried forward.
    frame = frame.loc[~frame.claim_id.isin(["SET-009", "SET-010"])].copy()
    frame["manuscript_location"] = frame.manuscript_location.astype(str).map(
        lambda value: f"V14 corresponding section (V13 ledger location: {value})"
    )
    frame["v14_evidence_update"] = "V13 audit result remains applicable; V14 did not change the numerical value."
    for claim_id, (value, note) in RESOLVED.items():
        mask = frame.claim_id.eq(claim_id)
        if mask.sum() != 1:
            raise RuntimeError(f"expected exactly one {claim_id} row")
        frame.loc[mask, "independently_recomputed_value"] = value
        frame.loc[mask, "difference"] = "consistent"
        frame.loc[mask, "rounding_status"] = "verified for V14"
        frame.loc[mask, "PASS_WARN_FAIL"] = "PASS"
        frame.loc[mask, "notes"] = note
        frame.loc[mask, "v14_evidence_update"] = note

    sim21 = frame.claim_id.eq("SIM-021")
    frame.loc[sim21, "independently_recomputed_value"] = (
        "Vacuum/air speed difference gives only +0.00626 deg at 20 deg and +0.01442 deg at 40 deg"
    )
    frame.loc[sim21, "difference"] = "+89,911 m/s (+0.0300%)"
    frame.loc[sim21, "rounding_status"] = "numerically negligible; wording remains slightly overbroad"
    frame.loc[sim21, "PASS_WARN_FAIL"] = "WARN"
    frame.loc[sim21, "notes"] = (
        "Say 'same planar-fit geometry and zenith cut' or harmonize the propagation-speed constant."
    )
    frame.loc[sim21, "v14_evidence_update"] = frame.loc[sim21, "notes"]

    additions = pd.DataFrame(
        [
            {
                "claim_id": "GEO-001",
                "manuscript_location": "V14 detector setup",
                "reported_value": "pentagonal footprint about 43 m^2",
                "quantity_meaning": "Shoelace area enclosed by the five measured detector centers",
                "authoritative_source_file": "V12 detector coordinates",
                "independently_recomputed_value": "42.946150395 m^2",
                "difference": "-0.053849605 m^2 from rounded 43 m^2",
                "rounding_status": "correctly rounded",
                "PASS_WARN_FAIL": "PASS",
                "notes": "Independent geometry calculation.",
                "v14_evidence_update": "New V14 claim verified in events_coordinates_exposure_geometry.json.",
            },
            {
                "claim_id": "GEO-002",
                "manuscript_location": "V14 detector setup",
                "reported_value": "five scintillators provide 1.25 m^2 active area",
                "quantity_meaning": "5 x 0.5 m x 0.5 m horizontal scintillator area",
                "authoritative_source_file": "Detector geometry description and GEANT4 metadata",
                "independently_recomputed_value": "1.25 m^2",
                "difference": "0",
                "rounding_status": "exact",
                "PASS_WARN_FAIL": "PASS",
                "notes": "Independent dimensional calculation.",
                "v14_evidence_update": "New V14 claim verified in events_coordinates_exposure_geometry.json.",
            },
            {
                "claim_id": "SIM-022",
                "manuscript_location": "V14 gamma-ray simulation, angular-resolution interpretation",
                "reported_value": "simulated r68 near 11 degrees is independently consistent with the empirical resolution",
                "quantity_meaning": "Independence of the simulated angular-resolution check",
                "authoritative_source_file": "pentagon_unthinned_effective_area.py and timing sensitivity audit",
                "independently_recomputed_value": (
                    "The undocumented default 1.85 ns station jitter alone gives r68=7.883, 8.153, and 7.715 deg "
                    "for ideal 0, 20, and 40 deg planar fronts"
                ),
                "difference": "The simulation imports a substantial assumed timing-smearing scale",
                "rounding_status": "numerical containment reproduces; independence wording is unsupported",
                "PASS_WARN_FAIL": "WARN",
                "notes": (
                    "Disclose and empirically justify 1.85 ns, show a jitter sensitivity envelope, and replace "
                    "'independently consistent' unless that value was fixed from an external calibration."
                ),
                "v14_evidence_update": "New issue found by the fresh timing-sensitivity audit.",
            },
        ]
    )
    frame = pd.concat([frame, additions], ignore_index=True)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUTPUT, index=False)
    counts = frame.PASS_WARN_FAIL.value_counts().to_dict()
    print(f"V14 claim groups: {len(frame)}; statuses: {counts}")
    print(frame.loc[frame.PASS_WARN_FAIL.ne("PASS"), ["claim_id", "PASS_WARN_FAIL", "notes"]].to_string(index=False))


if __name__ == "__main__":
    main()
