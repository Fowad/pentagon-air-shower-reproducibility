# Reproducibility scope

The package is intended to reproduce the numerical basis of the near-final V14 paper, not to perform new catalog fishing or extend the scientific question.

## Experimental branch

The authoritative event sample is rebuilt from the original 44 DAQ files. The V12 analysis code remains the primary reconstruction/map engine. A post-V12 archival layer then retains products that were calculated during the final scientific investigation but were not all written by the original V12 output routine: the 37-source first-LHAASO covariance/population tests, the focused 17-direction UHE family and 102 source-interval family, direct checkpointed LS I +61°303 phase permutations, the 8°/11° robustness ensemble, and live-exposure/fluence products.

## Simulation branch

The detector response is regenerated with geant4-pybind 0.1.3 / Geant4 11.4.1 and the documented 50x50x2 cm PVT + 1 mm iron geometry. The CORSIKA executable is the locked gamma-only application at commit `601fe3036725af876d6ab46a7cfe15f0c0987be5`, with SIBYLL 2.3d, UrQMD, SOPHIA and PROPOSAL. The science output workspace is separate from the historical workspace. Fixed seeds reproduce the 100/200/300 TeV, 0/20/40 degree response map and the limited azimuth study. Raw ground-particle Parquet files are retained.

## What “from scratch” means

Scientific data products are regenerated from their upstream physical inputs. Software installations do not need to be pointlessly rebuilt if their source/binary provenance passes the packaged audit. If CORSIKA is absent, the bootstrap can construct it from the locked source commit. This distinction avoids conflating software compilation with scientific independence.
