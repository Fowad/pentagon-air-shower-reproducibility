# Pentagon Air-Shower Directional Analysis and Reproducibility Code

Scientific analysis and reproducibility code for the submitted manuscript:

> **F. Motahari, S. Mortazavi Moghaddam, and M. Bahmanabadi, _Time-Resolved Directional Analysis and Search for Ultra-High-Energy Sources with a Pentagon Air Shower Array_ (submitted, 2026).**

This repository contains the code used to reconstruct and statistically analyze a year-spanning five-detector extensive-air-shower data set, together with the targeted CORSIKA/Geant4 response-study code and independent numerical-audit scripts. The public repository deliberately excludes raw/collaboration-controlled data and bulky generated shower products.

## What the code does

- parses the original variable-length detector stream from 44 observing runs;
- performs runwise timing calibration and five-detector direction reconstruction;
- applies the final zenith acceptance selection and celestial coordinate transforms;
- builds within-run time-randomized reference skies and independent calibration ensembles;
- carries out blind directional searches, exact source-centered tests, source-family tests, equal-live-time interval checks, and look-elsewhere calibration;
- archives covariance/permutation products, aperture-robustness checks, source-exposure calculations, and reproducibility manifests;
- models scintillator energy deposition with Geant4 and folds fixed-seed CORSIKA gamma showers through the five-detector geometry;
- provides independent scripts for event, coordinate, exposure, source-family, timing-jitter, and simulation cross-checks.

The definitive archival analysis used **532,527 selected events**, **44 runs**, and **5,032.673361 h** of live exposure, with **100,000 randomized reference skies** and **100,000 independent calibration skies**. These values are validation anchors, not hard-coded analysis inputs.

## Repository structure

```text
.
├── experimental/
│   ├── analysis/                 # V12 reconstruction/statistics + scientific audits
│   ├── frozen_inputs/            # small source lists + catalog provenance manifest
│   ├── reference_results/        # compact numerical comparison anchors
│   └── run_full_reproducibility.py
├── simulation/
│   ├── geant4_code/              # detector deposited-energy response
│   ├── support/                   # CORSIKA application + shower folding
│   ├── reference_results/        # protocol + compact reference summaries
│   └── *.command / *.py          # macOS/bootstrap and campaign launchers
├── audit/                         # independent final numerical-audit scripts
├── docs/                          # scope, numerical anchors, QA notes
├── scripts/repository_smoke_check.py
├── environment.yml
├── requirements.txt
├── DATA_AVAILABILITY.md
├── THIRD_PARTY_NOTICES.md
├── CONTRIBUTING.md
├── CITATION.cff
└── LICENSE
```

## Quick start: experimental branch

Create an environment:

```bash
conda env create -f environment.yml
conda activate pentagon-air-shower-repro
```

Then provide the 44 raw DAQ text files locally and run:

```bash
python experimental/analysis/Pentagon_Array_Directional_Analysis_V12.py \
  --data-dir /path/to/Data_Full_5_Center \
  --output-dir results/v12 \
  --reference-realizations 100000 \
  --look-elsewhere-realizations 100000 \
  --targeted-source-realizations 10000
```

The full scientific run is intentionally expensive. Reduced Monte Carlo counts are suitable only for development/smoke testing.

## Reproducibility and data policy

The public repository contains **no raw DAQ event stream** and **no raw CORSIKA shower files**. It also omits large generated checkpoint arrays and full byte copies of third-party catalogs. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

For exact catalog-level archival reproduction, the internal run used frozen 2026 catalog snapshots identified by SHA-256 in `experimental/frozen_inputs/catalog_snapshot_manifest_V12.csv`. Without that cache, V12 can query the official services and records the hashes of the new downloads.

## Static QA

Run:

```bash
python scripts/repository_smoke_check.py
```

This compiles all Python source files, parses JSON/YAML metadata, checks CSV readability, verifies that prohibited raw-data file types are absent, and syntax-checks shell launchers when `bash` is available. GitHub Actions runs the same check automatically.

## Simulation branch

The simulation branch requires additional external software and is documented separately in [simulation/README.md](simulation/README.md). It uses a threshold-response envelope rather than assigning event-by-event primary energies to the archived data.

## Citation

If you use the repository, cite the software metadata in `CITATION.cff` and the associated manuscript above. The repository's MIT license applies to project-authored code/documentation only; external software and catalog data retain their own terms.
