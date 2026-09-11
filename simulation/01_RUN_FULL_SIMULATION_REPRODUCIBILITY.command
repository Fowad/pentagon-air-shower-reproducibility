#!/bin/bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PENTAGON_CORSIKA_SOFTWARE_WORK="${PENTAGON_CORSIKA_SOFTWARE_WORK:-$HOME/Desktop/Pentagon_CORSIKA_Response_Stage1}"
export PENTAGON_REPRO_SIM_WORK="${PENTAGON_REPRO_SIM_WORK:-$HOME/Desktop/Pentagon_CORSIKA_FINAL_REPRODUCIBILITY}"
SOFT="$PENTAGON_CORSIKA_SOFTWARE_WORK"
WORK="$PENTAGON_REPRO_SIM_WORK"
if command -v caffeinate >/dev/null 2>&1; then caffeinate -dimsu -w $$ >/dev/null 2>&1 & fi
if [[ ! -x "$SOFT/corsika8-build/applications/c8_air_shower" ]]; then
  echo "[SETUP] Validated CORSIKA binary absent; bootstrapping first."
  bash "$ROOT/00_BOOTSTRAP_CORSIKA_IF_NEEDED.command"
fi
CONDA_EXE=""
for c in "$(command -v conda 2>/dev/null || true)" "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda"; do [[ -n "$c" && -x "$c" ]] && { CONDA_EXE="$c"; break; }; done
[[ -n "$CONDA_EXE" ]] || { echo "Conda not found after bootstrap"; exit 2; }
CONDA_BASE="$($CONDA_EXE info --base)"; source "$CONDA_BASE/etc/profile.d/conda.sh"; set +u; conda activate pentagon-corsika8; set -u
python - <<'PY' >/dev/null 2>&1 || python -m pip install --disable-pip-version-check "geant4-pybind==0.1.3"
import geant4_pybind
PY
python - <<'PY'
from importlib.metadata import version
import geant4_pybind
print('[GEANT4] geant4-pybind',version('geant4-pybind'))
PY
if [[ -e "$WORK" ]]; then
  echo "[SAFETY] $WORK already exists. This launcher will resume validated partial products rather than overwrite them."
else
  mkdir -p "$WORK"
fi
mkdir -p "$WORK/geant4_response"
LOG="$WORK/FULL_SIMULATION_REPRODUCIBILITY.log"
exec > >(tee -a "$LOG") 2>&1
# Regenerate detector response from first principles.
MIP="$WORK/geant4_response/mip_muon_10gev_vertical_steel_20k.npz"
if [[ ! -f "$MIP" ]]; then
  python "$ROOT/geant4_code/pentagon_scintillator_geant4.py" --particle mu- --energy-mev 10000 --angle-deg 0 --events 20000 --seed 420260820 --output "$MIP"
fi
python "$ROOT/geant4_code/derive_mip_calibration.py" --input "$MIP" --output "$WORK/geant4_response/mip_muon_10gev_vertical_steel_20k_summary.json"
python "$ROOT/geant4_code/run_geant4_response_grid.py" --driver "$ROOT/geant4_code/pentagon_scintillator_geant4.py" --output-dir "$WORK/geant4_response/particle_grid_v1" --events-em 3000 --events-other 2000 --workers 6
python "$ROOT/geant4_code/analyze_geant4_response_grid.py" --grid-dir "$WORK/geant4_response/particle_grid_v1"
python "$ROOT/geant4_code/run_geant4_high_energy_em_extension.py" --driver "$ROOT/geant4_code/pentagon_scintillator_geant4.py" --output-dir "$WORK/geant4_response/high_energy_em_extension_v1" --events 2000 --workers 1
# Run every CORSIKA shower and every response fold used in the supporting study.
python "$ROOT/full_gamma_campaign.py"
# Archive executing code and a hash manifest including raw particles.
rm -rf "$WORK/EXECUTING_CODE_SNAPSHOT"; cp -R "$ROOT" "$WORK/EXECUTING_CODE_SNAPSHOT"
python - "$WORK" <<'PY'
from pathlib import Path
import hashlib,json,sys,time
r=Path(sys.argv[1]); rows=[]
for p in sorted(r.rglob('*')):
 if p.is_file() and 'EXECUTING_CODE_SNAPSHOT' not in p.parts and p.name!='FINAL_SIMULATION_FILE_MANIFEST.json':
  h=hashlib.sha256();
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1<<20),b''): h.update(b)
  rows.append({'path':str(p.relative_to(r)),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
(r/'FINAL_SIMULATION_FILE_MANIFEST.json').write_text(json.dumps({'created_unix':time.time(),'files':rows},indent=2)+'\n')
print('[ARCHIVE] complete file manifest written:',len(rows),'files')
PY
echo "[COMPLETE] Full GEANT4 + CORSIKA supporting simulation reproduction finished: $WORK"
