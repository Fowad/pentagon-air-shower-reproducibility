# CORSIKA / Geant4 response study

This branch contains the targeted detector-response simulation supporting the manuscript. It is **not** an event-by-event energy reconstruction of the measured archive.

The simulation model includes a 50 x 50 x 2 cm PVT scintillator in a 1 mm iron enclosure, a Geant4 deposited-energy response model, and a fixed-seed CORSIKA 8 gamma-shower campaign folded across the measured five-detector geometry.

The detailed scientific protocol is in `reference_results/RESPONSE_STUDY_PROTOCOL.md`.

## Software

- Python 3.11
- `geant4-pybind==0.1.3` / Geant4 11.4.1 for the archived detector-response configuration
- CORSIKA 8, locked in the archival study to commit `601fe3036725af876d6ab46a7cfe15f0c0987be5`

CORSIKA and Geant4 are not vendored in this repository.

## Public-repository exclusions

Raw ground-particle Parquet files and large generated Geant4 response arrays are intentionally excluded. The code needed to regenerate them and small reference summaries/manifests are included.

The `.command` files were designed for the original macOS reproduction environment. For other systems, read the commands before running and adapt environment/bootstrap paths as needed.
