#!/bin/bash
set -Eeuo pipefail

# Pentagon paper: R18 runtime-resume + technical smoke test.
# R14/R15 proved that the paper-minimal CORSIKA executable builds and links.
# R16 then proved runtime dylib resolution and --help startup, but stopped on a false
# substring-based verification guard. This script deliberately DOES NOT rebuild CORSIKA. It validates the existing binary,
# repairs runtime search paths for the conda-forge Fortran dylibs + SIBYLL framework,
# then runs --help and one tiny 10 GeV gamma technical smoke shower.
# No scientific response inference is authorized here.

ROOT="$(cd "$(dirname "$0")" && pwd)"
WORK="${PENTAGON_CORSIKA_SOFTWARE_WORK:-$HOME/Desktop/Pentagon_CORSIKA_Response_Stage1}"
LOGDIR="$WORK/logs"
ENV_NAME="pentagon-corsika8"
CORSIKA_COMMIT="601fe3036725af876d6ab46a7cfe15f0c0987be5"
CUSTOM_CPP_SHA="3214d6e3886669e70788d864b780490c3af0fce05ab6f532ef26aad876604f08"
SRC="$WORK/corsika8-source"
BUILD="$WORK/corsika8-build"
BINARY="$BUILD/applications/c8_air_shower"
SMOKE="$WORK/technical_smoke_gamma_10GeV"
RUNTIME_ENV="$WORK/CORSIKA_RUNTIME_ENV.sh"
LOG="$LOGDIR/01_runtime_resume_and_smoke_R18.log"

mkdir -p "$LOGDIR"
exec > >(tee -a "$LOG") 2>&1

if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -dimsu -w $$ >/dev/null 2>&1 &
fi

make_diagnostic_zip() {
  local tag="${1:-diagnostic}"
  local out="$WORK/CORSIKA_STAGE1_R18_${tag}_DIAGNOSTIC.zip"
  rm -f "$out"
  python3 - "$WORK" "$out" <<'PY'
import sys,zipfile
from pathlib import Path
work=Path(sys.argv[1]); out=Path(sys.argv[2])
patterns=[
 'logs/*.log','logs/*.txt','BUILD_PROVENANCE.json','R18_RUNTIME_LINK_AUDIT.json','R18_REUSE_AUDIT.json','CORSIKA_RUNTIME_ENV.sh',
 'technical_smoke_gamma_10GeV/config.yaml','technical_smoke_gamma_10GeV/summary.yaml','technical_smoke_gamma_10GeV/primary/*.yaml',
 'corsika8-build/applications/CMakeFiles/c8_air_shower.dir/link.txt',
 'corsika8-source/applications/c8_air_shower.cpp','corsika8-source/corsika/modules/sibyll/Decay.hpp'
]
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for pattern in patterns:
        for p in work.glob(pattern):
            if p.is_file(): z.write(p,p.relative_to(work))
print(out)
PY
}

on_error() {
  code=$?
  echo
  echo "[FAILED] R18 runtime/smoke stage stopped with exit code $code."
  echo "Nothing from this failed stage is a scientific response result."
  make_diagnostic_zip "RUNTIME_OR_SMOKE_FAILURE" || true
  echo "Please send the R18 diagnostic ZIP from: $WORK"
  echo
  exit "$code"
}
trap on_error ERR

echo "============================================================"
echo " PENTAGON PAPER - R18 RUNTIME RESUME + SMOKE TEST"
echo "============================================================"
echo "Package: $ROOT"
echo "Workspace: $WORK"
echo "Locked CORSIKA commit: $CORSIKA_COMMIT"
echo

[[ "$(uname -s)" == "Darwin" ]] || { echo "This helper is for the user's macOS laptop."; exit 2; }
[[ "$(uname -m)" == "x86_64" ]] || { echo "R18 is reviewed only for the Intel x86_64 Mac."; exit 2; }
echo "macOS: $(sw_vers -productVersion 2>/dev/null || true)"
echo "Architecture: $(uname -m)"

# Locate and activate the existing isolated conda environment. No installs and no network.
CONDA_EXE=""
if command -v conda >/dev/null 2>&1; then
  CONDA_EXE="$(command -v conda)"
else
  for candidate in "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
    if [[ -x "$candidate" ]]; then CONDA_EXE="$candidate"; break; fi
  done
fi
[[ -n "$CONDA_EXE" ]] || { echo "Conda not found."; exit 2; }
CONDA_BASE="$($CONDA_EXE info --base)"
# shellcheck disable=SC1090
source "$CONDA_BASE/etc/profile.d/conda.sh"
set +u
conda activate "$ENV_NAME"
set -u
echo "Conda env: $CONDA_PREFIX"

# This is intentionally a resume script. Fail closed rather than rebuilding.
[[ -d "$SRC/.git" ]] || { echo "Existing CORSIKA source checkout missing. R18 will not rebuild it."; exit 3; }
[[ -x "$BINARY" ]] || { echo "Existing minimal CORSIKA binary missing. R18 will not rebuild it."; exit 3; }
[[ -f "$SRC/applications/c8_air_shower.cpp" ]] || { echo "Existing minimal application source missing."; exit 3; }
HEAD="$(git -C "$SRC" rev-parse HEAD)"
[[ "$HEAD" == "$CORSIKA_COMMIT" ]] || { echo "Source HEAD is not the locked CORSIKA commit: $HEAD"; exit 3; }
SRC_SHA="$(shasum -a 256 "$SRC/applications/c8_air_shower.cpp" | awk '{print $1}')"
PKG_SHA="$(shasum -a 256 "$ROOT/support/c8_pentagon_gamma_minimal.cpp" | awk '{print $1}')"
[[ "$PKG_SHA" == "$CUSTOM_CPP_SHA" ]] || { echo "Package minimal-source checksum mismatch."; exit 3; }
[[ "$SRC_SHA" == "$CUSTOM_CPP_SHA" ]] || { echo "Workspace application does not match the R14/R15 minimal source."; echo "workspace=$SRC_SHA"; echo "expected=$CUSTOM_CPP_SHA"; exit 3; }

LINKTXT="$BUILD/applications/CMakeFiles/c8_air_shower.dir/link.txt"
[[ -f "$LINKTXT" ]] || { echo "Generated final link.txt missing; cannot verify reused binary."; exit 3; }
# Require retained baseline components and reject broad model families in the actual final link command.
for token in SIBYLL23d Sophia UrQMD PROPOSAL; do
  grep -qi "$token" "$LINKTXT" || { echo "Reuse audit failed: baseline token absent from final link command: $token"; exit 3; }
done
if grep -Eqi 'QGSJet|EPOS_LHC|EPOS_LHCR|Tauola|Pythia|CONEX' "$LINKTXT"; then
  echo "Reuse audit failed: broad-model family appears in final link command."
  grep -Ei 'QGSJet|EPOS_LHC|EPOS_LHCR|Tauola|Pythia|CONEX' "$LINKTXT" || true
  exit 3
fi

BINARY_SHA_BEFORE="$(shasum -a 256 "$BINARY" | awk '{print $1}')"
python - "$WORK/R18_REUSE_AUDIT.json" "$HEAD" "$SRC_SHA" "$BINARY_SHA_BEFORE" "$LINKTXT" <<'PY'
from pathlib import Path
import json,sys,time
out=Path(sys.argv[1])
out.write_text(json.dumps({
 'corsika_commit':sys.argv[2],
 'minimal_source_sha256':sys.argv[3],
 'binary_sha256_before_runtime_rpath_repair':sys.argv[4],
 'link_txt':sys.argv[5],
 'broad_models_absent_from_final_link':True,
 'baseline_components_verified':['SIBYLL23d','Sophia','UrQMD','PROPOSAL'],
 'created_unix':time.time()
},indent=2,sort_keys=True)+'\n')
PY
echo "[REUSE] Existing paper-minimal binary/source provenance PASSED."
echo "[REUSE] No CORSIKA rebuild will be performed."
echo "[REUSE] Binary SHA-256 before runtime repair: $BINARY_SHA_BEFORE"

command -v gfortran >/dev/null 2>&1 || { echo "Active gfortran missing."; exit 4; }
command -v otool >/dev/null 2>&1 || { echo "otool missing."; exit 4; }
command -v install_name_tool >/dev/null 2>&1 || { echo "install_name_tool missing."; exit 4; }
FC_REAL="$(command -v gfortran)"

# Resolve runtime dylibs WITHOUT Bash arrays. macOS ships an older Bash where an
# empty array expanded under `set -u` can abort; R15 died here before any runtime test.
FORTRAN_DYLD_PATH=""
resolve_runtime_lib() {
  local lib="$1"
  local p=""
  local d=""
  p="$($FC_REAL -print-file-name="$lib" 2>/dev/null || true)"
  if [[ -z "$p" || "$p" == "$lib" || ! -e "$p" ]]; then
    p="$(find "$CONDA_PREFIX" -type f -name "$lib" -print -quit 2>/dev/null || true)"
  fi
  [[ -n "$p" && -e "$p" ]] || { echo "Unable to locate required Fortran runtime dylib: $lib"; return 1; }
  d="$(cd "$(dirname "$p")" && pwd)"
  case ":$FORTRAN_DYLD_PATH:" in
    *":$d:"*) ;;
    *) FORTRAN_DYLD_PATH="${FORTRAN_DYLD_PATH:+$FORTRAN_DYLD_PATH:}$d" ;;
  esac
  echo "[RUNTIME] $lib -> $p"
}

resolve_runtime_lib libgfortran.5.dylib
resolve_runtime_lib libgcc_s.1.1.dylib
resolve_runtime_lib libquadmath.0.dylib
[[ -n "$FORTRAN_DYLD_PATH" ]] || { echo "Fortran runtime path discovery produced no directories."; exit 4; }

DATA_DIR="$SRC/modules/data"
SIBYLL_DIR="$BUILD/modules/sibyll"
[[ -f "$DATA_DIR/GeoMag/IGRF13.COF" ]] || { echo "Missing CORSIKA geomagnetic data file."; exit 4; }
[[ -d "$SIBYLL_DIR/SIBYLL23d.framework" ]] || { echo "Built SIBYLL framework missing."; exit 4; }

python - "$RUNTIME_ENV" "$DATA_DIR" "$SIBYLL_DIR" "$FORTRAN_DYLD_PATH" <<'PY'
from pathlib import Path
import shlex,sys
out=Path(sys.argv[1]); data=Path(sys.argv[2]).resolve(); sib=Path(sys.argv[3]).resolve(); fdirs=sys.argv[4]
lines=[
 '#!/bin/bash',
 '# Generated by R18 runtime-resume. Runtime paths only; no scientific response settings.',
 'export CORSIKA_DATA='+shlex.quote(str(data)),
 'export DYLD_FRAMEWORK_PATH='+shlex.quote(str(sib))+'${DYLD_FRAMEWORK_PATH:+:$DYLD_FRAMEWORK_PATH}',
 'export DYLD_LIBRARY_PATH='+shlex.quote(fdirs)+'${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}',
 'export DYLD_FALLBACK_LIBRARY_PATH='+shlex.quote(fdirs)+'${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}',
]
out.write_text('\n'.join(lines)+'\n')
PY
chmod +x "$RUNTIME_ENV"
# shellcheck disable=SC1090
source "$RUNTIME_ENV"
echo "[RUNTIME] CORSIKA_DATA=$CORSIKA_DATA"
echo "[RUNTIME] DYLD_FRAMEWORK_PATH=$DYLD_FRAMEWORK_PATH"
echo "[RUNTIME] DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH"

# Add durable LC_RPATH entries to the already-built binary only when truly absent.
# IMPORTANT: do not use an `otool | ... | grep -q` pipeline here. With `set -o pipefail`,
# grep -q can exit after a match and give an upstream process SIGPIPE, falsely making the
# whole pipeline fail. R18 hit exactly that false-negative and attempted a duplicate rpath.
RPATHS="$SIBYLL_DIR:$FORTRAN_DYLD_PATH"
OLDIFS="$IFS"
IFS=':'
for d in $RPATHS; do
  [[ -n "$d" ]] || continue
  if python - "$BINARY" "$d" <<'PYRPATH'
import subprocess,sys
binary,want=sys.argv[1],sys.argv[2]
text=subprocess.check_output(['otool','-l',binary], text=True, stderr=subprocess.STDOUT)
paths=[]
lines=text.splitlines()
for i,line in enumerate(lines):
    if line.strip() == 'cmd LC_RPATH':
        # The LC_RPATH load command contains a later line like:
        #     path /some/directory (offset 12)
        for j in range(i+1, min(i+6, len(lines))):
            t=lines[j].strip()
            if t.startswith('path '):
                val=t[5:]
                if ' (offset ' in val:
                    val=val.split(' (offset ',1)[0]
                paths.append(val)
                break
if want in paths:
    print(f'[RUNTIME] LC_RPATH already present: {want}')
    raise SystemExit(0)
raise SystemExit(1)
PYRPATH
  then
    :
  else
    echo "[RUNTIME] Adding LC_RPATH: $d"
    install_name_tool -add_rpath "$d" "$BINARY"
    # Verify the add really took effect; fail closed if it did not.
    python - "$BINARY" "$d" <<'PYRPATHVERIFY'
import subprocess,sys
binary,want=sys.argv[1],sys.argv[2]
text=subprocess.check_output(['otool','-l',binary], text=True, stderr=subprocess.STDOUT)
paths=[]
lines=text.splitlines()
for i,line in enumerate(lines):
    if line.strip() == 'cmd LC_RPATH':
        for j in range(i+1, min(i+6, len(lines))):
            t=lines[j].strip()
            if t.startswith('path '):
                val=t[5:]
                if ' (offset ' in val:
                    val=val.split(' (offset ',1)[0]
                paths.append(val)
                break
if want not in paths:
    raise SystemExit(f'LC_RPATH verification failed after add: {want}')
print(f'[RUNTIME] LC_RPATH add verified: {want}')
PYRPATHVERIFY
  fi
done
IFS="$OLDIFS"

# If the local linker produced a code signature, refresh it after Mach-O modification.
if command -v codesign >/dev/null 2>&1 && codesign -dv "$BINARY" >/dev/null 2>&1; then
  echo "[RUNTIME] Refreshing ad-hoc code signature after LC_RPATH update."
  codesign --force --sign - "$BINARY" >/dev/null 2>&1
fi

otool -L "$BINARY" | tee "$LOGDIR/c8_air_shower_otool_L_R18.txt"
otool -l "$BINARY" > "$LOGDIR/c8_air_shower_otool_l_R18.txt"

python - "$LOGDIR/c8_air_shower_otool_L_R18.txt" "$FORTRAN_DYLD_PATH" "$WORK/R18_RUNTIME_LINK_AUDIT.json" <<'PY'
from pathlib import Path
import json,sys,time
otool=Path(sys.argv[1]).read_text(errors='replace')
dirs=[Path(x) for x in sys.argv[2].split(':') if x]
required=['libgfortran.5.dylib','libgcc_s.1.1.dylib','libquadmath.0.dylib']
report={'runtime_dirs':[str(x) for x in dirs],'required':{},'otool_dependencies':otool.splitlines()[1:],'created_unix':time.time()}
for lib in required:
    hits=[str(d/lib) for d in dirs if (d/lib).exists()]
    report['required'][lib]={'resolved_candidates':hits,'appears_in_otool':any(lib in line for line in otool.splitlines())}
    if not hits:
        Path(sys.argv[3]).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
        raise SystemExit('Required Fortran runtime dylib not resolvable: '+lib)
Path(sys.argv[3]).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print('[RUNTIME] Fortran dylib resolution audit PASSED.')
PY

BINARY_SHA_AFTER="$(shasum -a 256 "$BINARY" | awk '{print $1}')"
echo "[RUNTIME] Binary SHA-256 after LC_RPATH repair: $BINARY_SHA_AFTER"

# Startup probe first. On any failure, capture dyld search and loaded-library diagnostics.
set +e
"$BINARY" --help > "$LOGDIR/c8_air_shower_help_R18.txt" 2>&1
HELP_RC=$?
set -e
if (( HELP_RC != 0 )); then
  echo "[RUNTIME] --help failed with exit code $HELP_RC; capturing dyld diagnostics."
  set +e
  DYLD_PRINT_SEARCHING=1 DYLD_PRINT_LIBRARIES=1 "$BINARY" --help > "$LOGDIR/c8_air_shower_dyld_R18.txt" 2>&1
  set -e
  echo "--- --help output tail ---"
  tail -n 100 "$LOGDIR/c8_air_shower_help_R18.txt" 2>/dev/null || true
  echo "--- dyld diagnostic tail ---"
  tail -n 160 "$LOGDIR/c8_air_shower_dyld_R18.txt" 2>/dev/null || true
  false
fi

echo "[RUNTIME] --help startup probe PASSED."
python - "$LOGDIR/c8_air_shower_help_R18.txt" <<'PYHELP'
from pathlib import Path
import re,sys
text=Path(sys.argv[1]).read_text(errors='replace')
# Parse exact long-option tokens. R16 incorrectly used substring grep, so the
# legitimate --hadronModelTransitionEnergy matched the forbidden --hadronModel.
opts=set(re.findall(r'--[A-Za-z][A-Za-z0-9-]*', text))
required={
    '--energy','--geomagnetic-year','--site-latitude','--site-longitude',
    '--emcut','--hadcut','--mucut','--hadronModelTransitionEnergy'
}
missing=sorted(required-opts)
if missing:
    raise SystemExit('Required minimal option(s) absent: '+', '.join(missing))
forbidden={'--hadronModel','--pdg'}
present=sorted(forbidden & opts)
if present:
    raise SystemExit('Forbidden broad selector(s) present: '+', '.join(present))
print('[RUNTIME] Exact-option minimal-app verification PASSED.')
print('[RUNTIME] Fixed gamma primary; fixed SIBYLL/UrQMD model chain confirmed by source/link provenance.')
PYHELP

# Technical smoke only. Not a paper result.
rm -rf "$SMOKE"
echo
echo "[SMOKE] Running one 10 GeV gamma technical shower with the minimal baseline chain."
set +e
"$BINARY" \
  --energy 10 --zenith 20 --azimuth 0 \
  --observation-level 1200 \
  --geomagnetic-year 2017.0 --site-latitude 35.704 --site-longitude 51.351 \
  --emcut 0.0005 --hadcut 0.3 --mucut 0.3 \
  --emthin 1e-6 --max-weight 1 --nevent 1 --seed 2026082001 \
  --verbosity info --filename "$SMOKE" > "$LOGDIR/technical_smoke_corsika_R18.log" 2>&1
SMOKE_RC=$?
set -e
if (( SMOKE_RC != 0 )); then
  echo "[SMOKE] Technical shower failed with exit code $SMOKE_RC."
  tail -n 160 "$LOGDIR/technical_smoke_corsika_R18.log" 2>/dev/null || true
  false
fi

python - "$SMOKE" <<'PY'
from pathlib import Path
import sys,yaml
import pyarrow.parquet as pq
p=Path(sys.argv[1])
if not (p/'summary.yaml').is_file(): raise SystemExit('smoke summary.yaml missing')
summary=yaml.safe_load((p/'summary.yaml').read_text())
if int(summary.get('showers',-1)) != 1: raise SystemExit(f'unexpected smoke shower count: {summary}')
q=p/'particles'/'particles.parquet'
if not q.is_file(): raise SystemExit('smoke ground-particle parquet missing')
meta=pq.read_metadata(q)
print(f'[SMOKE] Valid minimal CORSIKA output: showers=1, ground_rows={meta.num_rows}')
PY

python - "$WORK" "$SRC" "$BINARY" "$ROOT" "$BINARY_SHA_BEFORE" "$BINARY_SHA_AFTER" <<'PY'
from pathlib import Path
import hashlib,json,platform,subprocess,sys,time
work,src,binary,root=map(Path,sys.argv[1:5]); before=sys.argv[5]; after=sys.argv[6]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def cmd(*a):
    try:return subprocess.check_output(a,text=True,stderr=subprocess.STDOUT).strip()
    except Exception as e:return f'ERROR: {e}'
prov={
 'protocol':'PENTAGON-CORSIKA-PAPER-MINIMAL-STAGE1-R18-RUNTIME-RESUME',
 'purpose':'runtime repair + technical smoke only; no scientific response inference',
 'corsika_commit':cmd('git','-C',str(src),'rev-parse','HEAD'),
 'binary':str(binary),
 'binary_sha256_before_runtime_rpath_repair':before,
 'binary_sha256_after_runtime_rpath_repair':after,
 'minimal_source_sha256':sha(root/'support'/'c8_pentagon_gamma_minimal.cpp'),
 'retained_physics':['gamma primary','SIBYLL-2.3d high-energy hadronic','UrQMD low-energy hadronic','SIBYLL decays','PROPOSAL electromagnetic transport','SOPHIA resonance photonuclear','USStdBK atmosphere','IGRF-13','ground-particle Parquet output'],
 'stage02_authorized':False,
 'created_unix':time.time(),'platform':platform.platform(),
}
(work/'BUILD_PROVENANCE.json').write_text(json.dumps(prov,indent=2,sort_keys=True)+'\n')
PY

echo
echo "============================================================"
echo " BUILD + SMOKE TEST PASSED"
echo "============================================================"
echo "R18 reused the already-built PAPER-MINIMAL executable; no CORSIKA rebuild occurred."
echo "No scientific response has been inferred from the technical smoke shower."
echo "DO NOT run Stage 02."
echo "Next action: return to ChatGPT for the fresh manuscript-driven Stage 02 audit."
echo "Binary: $BINARY"
echo "Binary SHA-256 after runtime repair: $BINARY_SHA_AFTER"
echo
