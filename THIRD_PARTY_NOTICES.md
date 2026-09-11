# Third-party software and data

The MIT license in this repository applies only to code and documentation authored for this project. It does **not** relicense third-party software or external scientific databases.

## External software

The simulation branch interfaces with **CORSIKA 8**, **Geant4** (through `geant4-pybind`), and their own dependencies. Those projects are not vendored here and remain subject to their upstream licenses and citation requirements. The scientific Python stack (`NumPy`, `SciPy`, `pandas`, `Matplotlib`, `Astropy`, `healpy`, `PyYAML`, `PyArrow`) likewise retains its own licensing.

## External catalogs

The analysis queries public catalogs from services including HEASARC, Fermi, VizieR/CDS, H.E.S.S., HAWC, and the McGill magnetar catalog. The public repository includes only the **snapshot manifest** (source URLs, hashes, and row counts) used to document provenance. Full frozen third-party catalog byte snapshots from the internal archival package are intentionally not redistributed here.

## CORSIKA application source

`simulation/support/c8_pentagon_gamma_minimal.cpp` is a project-specific application built against the CORSIKA 8 API. Users must obtain and build CORSIKA 8 separately and comply with its upstream terms.
