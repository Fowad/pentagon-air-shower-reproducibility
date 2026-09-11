# Pentagon-array CORSIKA/GEANT4 response-study protocol

Locked before inspection of any source-yield or source-significance result.  Pilot
showers may be used to choose computational parameters, but not to tune the
detector response toward a desired energy scale or source conclusion.

## Scientific scope

The study asks what five-fold response is physically plausible for primary gamma
rays from 0.03 to 10 PeV.  It is not an energy reconstruction of the archived
events.  The discriminator settings, PMT gains, and per-station calibration are
not preserved, so no single absolute trigger threshold can be inferred from the
TAC-only archive.  All quantitative results must therefore remain an explicit
threshold-response envelope.

The response study is independent of the catalog/source tests already performed.
No simulated response parameter may be selected because it makes a catalog source
look more or less significant.

## Fixed detector and site model

- Five identical 50 x 50 x 2 cm PVT scintillators in 1 mm iron enclosures.
- Detector coordinates (m): (-2.7586,-3.6756), (-4.3133,1.0717),
  (-0.2780,4.0196), (3.7707,1.0942), (2.2376,-3.6617).
- Observation level: 1200 m above sea level.
- Site: 35.704 deg N, 51.351 deg E.
- Geomagnetic field: IGRF-13 evaluated at decimal year 2017.0.  In the
  north-west-up axes used by the array this is
  (27.9223, -2.2869, -39.3339) microtesla at 1200 m.
- Atmosphere: CORSIKA USStdBK.  The lack of a measured Tehran atmospheric
  profile is a named response systematic, not a fitted parameter.
- Hardware coincidence: all five stations within 200 ns.  The archived
  theta < 45 deg selection is applied after any simulated direction audit.

## Software identity

- CORSIKA 8 source commit: 601fe3036725af876d6ab46a7cfe15f0c0987be5.
- CORSIKA release string printed by the executable: 0.1.99.
- High-energy interaction baseline: SIBYLL 2.3d; QGSJet-II.04 is the
  predeclared model cross-check.
- Low-energy hadronic model: UrQMD 1.3.1.
- GEANT4 version: 11.4.1 through geant4-pybind 0.1.3; FTFP_BERT detector
  transport with a 0.1 mm production cut.
- The locally modified CORSIKA application replaces the stock example's
  generic 50 microtesla horizontal field with the recorded IGRF-13 site field.

Exact executable, source, command-line, random-seed, and output hashes are to be
written to the final manifest.

## Detector-response envelope

The GEANT4 vertical 10 GeV muon calibration fixes only a simulation unit:
1 MIP = 3.5785 MeV, the fitted most-probable deposited energy in the active
plastic.  It is not claimed to be the historical discriminator threshold.

Five-fold response is evaluated at deposited-energy thresholds of 0.25, 0.5,
1.0, and 2.0 MIP per station.  This range is retained in full; no single curve
will be selected after seeing an astrophysical result.  Optical collection,
PMT gain, and discriminator uncertainty are represented by this envelope.

GEANT4 ground-particle response tables will include photons, electrons,
positrons, muons, protons, neutrons, and charged pions over the kinetic-energy
and incidence-angle support found in the CORSIKA pilots.  Interpolation is in
log kinetic energy and secant angle.  Extrapolation outside the simulated table
must be flagged and bounded by an endpoint response.

## Shower grid

Primary-gamma grid:

- energies: 0.03, 0.10, 0.30, 1.0, 3.0, and 10 PeV;
- zenith angles: 0, 20, and 40 deg;
- azimuth: 0 deg for the main grid, with a 90 deg cross-check at 20 and 40 deg;
- independent seeds and first-interaction depths retained per shower.

Proton and helium simulations are response/rate controls, not gamma/hadron
classification.  They use the same energy and zenith support where computing
precision permits, with at least 0.10, 0.30, 1.0, and 3.0 PeV at 20 deg.

## Thinning and statistical stopping

Unthinned showers are the reference wherever computationally feasible.  Thinned
particles are never interpreted as simultaneous co-located physical particles.
Any accelerated calculation must convert weighted output into a local particle
density/compound-Poisson response and validate it against unthinned showers.

At 0.30 PeV and 20 deg, collect at least 12 unthinned gamma showers and continue,
up to 60, until the bootstrap 68% interval half-width of the median effective
area is no more than 25% for the 0.5- and 1-MIP curves, provided those curves are
nonzero.  Other unthinned grid points use a minimum of eight and maximum of 30
showers.  An accelerated response curve is quantitatively admissible only if its
median effective area differs from the unthinned reference by no more than 20%
at both 0.5 and 1 MIP and its shower-to-shower spread is reproduced to within
30%.  Failure of this check restricts thinning to qualitative diagnostics.

## Core translation and five-fold response

For every shower, translate the measured five-detector geometry across the
ground-particle footprint.  A ground particle intersects a detector only when
its trajectory crosses that detector's 0.5 x 0.5 m footprint.  Its deposited
energy is sampled from the applicable GEANT4 response table.  Station deposits
are formed with a 10 ns pulse-overlap window, with 5 ns and 20 ns retained as
predeclared electronics-shaping cross-checks.  The four fixed threshold
definitions and the 200 ns five-fold coincidence are then applied.

The core-integration grid is refined until changing its spacing changes the
effective area by less than 5%.  The integration boundary is enlarged until no
accepted core lies on its outer two grid cells.  Response randomization is
repeated independently so that GEANT4 sampling error is separated from genuine
shower-to-shower variation.

## Validation and permissible conclusions

Required checks are:

1. energy accounting, particle weights, cuts, and ground-file completion;
2. GEANT4 MIP deposit agreement with the expected approximately 4 MeV mean loss
   through 2 cm of plastic;
3. convergence in core-grid spacing and boundary;
4. unthinned-versus-accelerated response agreement under the stopping rule;
5. stability under SIBYLL 2.3d versus QGSJet-II.04, azimuth, and response seed;
6. qualitative consistency of the simulated zenith trend with the archived
   theta distribution, without fitting an unknown threshold to force agreement.

A manuscript figure is warranted only if the response envelope supports a stable
scientific statement across the declared thresholds and validation variants.
Otherwise the simulation remains a technical audit and cannot be used to assign
a calibrated energy threshold, effective area, or source-flux limit to the data.

## Pre-result clarification: site coordinates and transport cuts

This clarification was logged while the 0.30 PeV unthinned reference showers
were still running and before any unthinned effective area or source yield was
available.

The definitive V12 celestial conversion uses 35 deg 43 min N, 51 deg 20 min E.
The first reference-shower command used the decimal laboratory coordinates
35.704 deg N, 51.351 deg E for its IGRF evaluation.  The two locations are only
about 2 km apart.  Source-transit exposure will use the exact V12 constants.  A
site-field calculation at both coordinate pairs is required; the reference
shower need not be repeated only if the component-wise field change is
negligible relative to the other response systematics.

The baseline CORSIKA transport cuts are its documented application defaults:
0.5 MeV for photons/electrons/positrons and 0.3 GeV for hadrons and muons.  These
cuts are to be recorded explicitly.  Because several individually subthreshold
particles can sum inside a station, at least one lower-cut 0.30 PeV shower must
be compared with the baseline.  The lower-cut check uses 0.05 MeV for the
electromagnetic component and 0.02 GeV for hadrons and muons.  The GEANT4 table
must be extended below 0.5 MeV before those particles are folded; clipping them
to the 0.5 MeV endpoint is not permitted.  The cut comparison is a detector-
response systematic and may not be selected according to a source result.
