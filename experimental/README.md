# Experimental analysis

The primary analysis entry point is:

```text
analysis/Pentagon_Array_Directional_Analysis_V12.py
```

It parses the variable-length five-detector DAQ stream, applies runwise timing calibration, reconstructs directions with a plane fit, applies the zenith selection, transforms event directions to celestial coordinates, builds time-randomized reference skies, calibrates look-elsewhere effects, and performs directed source-family tests.

## Full scientific settings

```bash
python experimental/analysis/Pentagon_Array_Directional_Analysis_V12.py \
  --data-dir /path/to/Data_Full_5_Center \
  --output-dir results/v12 \
  --reference-realizations 100000 \
  --look-elsewhere-realizations 100000 \
  --targeted-source-realizations 10000
```

The 44 raw DAQ files are required and are not in this repository.

## Development/smoke settings

For code-path testing only (not scientific reproduction), the realization counts can be reduced substantially. Do not report probabilities from a reduced-trial smoke run as manuscript results.

## Catalog provenance

`frozen_inputs/catalog_snapshot_manifest_V12.csv` records the official sources and SHA-256 identities of the frozen public catalog snapshots used by the archival run. The byte snapshots themselves are not redistributed in the public repository. When no frozen cache is preseeded in the output directory, V12 downloads the current official catalog versions and records new hashes.

## Orchestrated driver

`run_full_reproducibility.py` provides explicit `--data-dir`, `--output-dir`, and optional `--catalog-cache-zip` arguments while preserving the core V12 entry point.
