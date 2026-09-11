# Data availability and public-repository scope

This repository is intentionally **code-first** and contains no raw experimental event stream.

## Not included

- the 44 original DAQ text files;
- reconstructed per-event archives derived from those files;
- raw CORSIKA ground-particle Parquet files;
- large Geant4/CORSIKA Monte Carlo checkpoint arrays (`.npz`, `.parquet`);
- full frozen byte copies of third-party astronomical catalogs.

These exclusions keep collaboration-controlled/raw data and bulky generated products out of a public code repository.

## Included

- the event-reconstruction and directional-analysis source code;
- downstream source-family, permutation, aperture-robustness, exposure, and audit scripts;
- CORSIKA/Geant4 response-study source code and launchers;
- small derived source lists and reference JSON/CSV summaries;
- a catalog snapshot manifest recording source URLs, hashes, and row counts;
- numerical anchors and reproducibility documentation.

## Reproducing the experimental branch

Provide a directory containing the original 44 DAQ text files and run the analysis with `--data-dir`. For an archival byte-identical catalog preflight, provide the separately retained `catalog_cache_V12.zip`; otherwise the code can query the official catalog services and records the newly downloaded snapshot hashes.

Access to collaboration-controlled raw data, if appropriate, is handled separately from this public repository and is subject to the authors' and host laboratory's data-sharing policy.
