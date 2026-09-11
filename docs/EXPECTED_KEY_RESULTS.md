# Expected scientific anchors for the reproduction

These are comparison targets, not inputs to the calculation. The new run should
produce them from upstream data/simulations, subject to the exact stored seeds
and the printed rounding used by the manuscript.

## Experimental archive
- Accepted reconstructed events: 532,527
- DAQ runs: 44
- Live exposure: 5032.673361 h
- Accepted-event archive identity expected from the authoritative reconstruction:
  `f7f370ad96d69ef9a793472490b4f10f4deff4e6c62cd2211527f1c65e245a6b`
- Reference skies: 100,000
- Independent calibration skies: 100,000

## Blind directional search
- Full-archive high-occupancy maximum: Z about 2.565; map-wide p about 0.7442
- Strongest equal-live interval maximum: Z about 3.938; interval map p about 0.01568
- Six-interval family calibration: p about 0.09012

## Exact source families
The reproduction branch is specifically designed to retain the complete null
products behind the manuscript values, rather than only the final probabilities.
Expected manuscript-level anchors include:
- 37 first-LHAASO >100-TeV directions: family p about 0.718
- Equal-source combined statistic: p about 0.674
- 100-TeV-flux-weighted combined statistic: p about 0.705
- Focused 17-direction >=~0.3-PeV family: full family p about 0.356
- 17 x 6 source-interval family: p about 0.138

## Gamma-response simulation
At 20 degrees and nominal 0.5-MIP / 10-ns response assumptions:
- 100 TeV: A_resp roughly 0.75-6.5 m^2 (very few accepted core positions)
- 200 TeV: A_resp roughly 170-763 m^2
- 300 TeV: A_resp roughly 345-1250 m^2
- Tested 40-degree showers through 300 TeV: no accepted fivefold response
- Substantially populated 200-300-TeV configurations: r68 roughly 11 degrees

The simulation reproduction should also retain raw shower Parquet files, the
GEANT4 MIP calibration and response-grid manifests, all fold replicas, and the
complete file-hash manifest.
