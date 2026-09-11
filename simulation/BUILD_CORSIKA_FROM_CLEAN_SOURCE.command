#!/bin/bash
set -Eeuo pipefail

# Pentagon paper: minimal local CORSIKA 8 build + technical smoke test.
# This stage deliberately builds only the baseline physics chain needed for the manuscript:
# gamma + SIBYLL 2.3d + UrQMD + SOPHIA + PROPOSAL + ground-particle output.
# It produces NO scientific response result; Stage 02 is intentionally not authorized yet.

ROOT="$(cd "$(dirname "$0")" && pwd)"
WORK="${PENTAGON_CORSIKA_SOFTWARE_WORK:-$HOME/Desktop/Pentagon_CORSIKA_Response_Stage1}"
LOGDIR="$WORK/logs"
ENV_NAME="pentagon-corsika8"
CORSIKA_URL="https://gitlab.iap.kit.edu/AirShowerPhysics/corsika.git"
CORSIKA_COMMIT="601fe3036725af876d6ab46a7cfe15f0c0987be5"
CUSTOM_CPP_SHA="3214d6e3886669e70788d864b780490c3af0fce05ab6f532ef26aad876604f08"
SRC="$WORK/corsika8-source"
BUILD="$WORK/corsika8-build"
BINARY="$BUILD/applications/c8_air_shower"
SMOKE="$WORK/technical_smoke_gamma_10GeV"
LOG="$LOGDIR/01_build_and_smoke.log"

mkdir -p "$LOGDIR"
exec > >(tee -a "$LOG") 2>&1

if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -dimsu -w $$ >/dev/null 2>&1 &
fi

make_diagnostic_zip() {
  local tag="${1:-diagnostic}"
  local out="$WORK/CORSIKA_STAGE1_${tag}_DIAGNOSTIC.zip"
  rm -f "$out"
  python3 - "$WORK" "$out" <<'PY'
import os, sys, zipfile
from pathlib import Path
work=Path(sys.argv[1]); out=Path(sys.argv[2])
patterns=[
    'logs/*.log','logs/c8_air_shower_otool_L.txt','BUILD_PROVENANCE.json','R14_SIBYLL_DECAY_DEFAULT_AUDIT.json','R14_MODULE_PRUNE.json','R14_SPDLOG_PORTABILITY_PATCH.json','R14_SIBYLL_LINK_PATCH.json','CORSIKA_RUNTIME_ENV.sh',
    'technical_smoke_gamma_10GeV/config.yaml','technical_smoke_gamma_10GeV/summary.yaml',
    'technical_smoke_gamma_10GeV/primary/*.yaml',
    'corsika8-build/CMakeCache.txt','corsika8-build/**/link.txt','corsika8-build/**/flags.make',
    'corsika8-source/applications/CMakeLists.txt','corsika8-source/applications/c8_air_shower.cpp','corsika8-source/corsika/modules/sibyll/Decay.hpp',
    'corsika8-source/modules/sibyll/CMakeLists.txt','corsika8-source/modules/sophia/CMakeLists.txt','corsika8-source/modules/urqmd/CMakeLists.txt',
    'corsika8-source/modules/pythia8/CMakeLists.txt','corsika8-source/modules/tauola/CMakeLists.txt','corsika8-source/modules/conex/CMakeLists.txt',
    'corsika8-source/modules/epos-lhcr/CMakeLists.txt','corsika8-source/modules/epos-lhc/CMakeLists.txt','corsika8-source/modules/qgsjetII/CMakeLists.txt','corsika8-source/modules/qgsjetIII/CMakeLists.txt',
    'corsika8-source/corsika/detail/framework/core/SpdlogSpecializations.inl','macos_compile_probes/*','application_preflight/*'
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
  echo "[FAILED] Build/smoke stage stopped with exit code $code."
  echo "Nothing from this failed stage is a scientific response result."
  make_diagnostic_zip "BUILD_OR_SMOKE_FAILURE" || true
  echo "Please send the diagnostic ZIP from: $WORK"
  echo
  exit "$code"
}
trap on_error ERR

echo "============================================================"
echo " PENTAGON PAPER - MINIMAL CORSIKA 8 BUILD + SMOKE TEST"
echo "============================================================"
echo "Package: $ROOT"
echo "Workspace: $WORK"
echo "Locked CORSIKA commit: $CORSIKA_COMMIT"
echo

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This helper was prepared for the user's macOS laptop."
  echo "Detected: $(uname -s)"
  exit 2
fi

ARCH="$(uname -m)"
echo "macOS: $(sw_vers -productVersion 2>/dev/null || true)"
echo "Architecture: $ARCH"
if [[ "$ARCH" != "x86_64" ]]; then
  echo "WARNING: this package was designed and reviewed for the Intel MacBook (x86_64)."
  echo "The upstream CORSIKA project has documented separate macOS/arm64 build issues."
  exit 2
fi

# Require conservative free disk for source/build. No production-shower allocation is authorized here.
FREE_KB=$(df -Pk "$HOME/Desktop" | awk 'NR==2 {print $4}')
FREE_GB=$((FREE_KB/1024/1024))
echo "Free Desktop disk space: approximately ${FREE_GB} GiB"
if (( FREE_GB < 20 )); then
  echo "At least 20 GiB free is required even to enter the disk-recovery/build preflight."
  exit 2
fi

# Locate the user's conda installation without modifying the existing eas environment.
CONDA_EXE=""
if command -v conda >/dev/null 2>&1; then
  CONDA_EXE="$(command -v conda)"
else
  for candidate in \
    "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" \
    "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
    "/opt/anaconda3/bin/conda" "/usr/local/anaconda3/bin/conda"; do
    if [[ -x "$candidate" ]]; then CONDA_EXE="$candidate"; break; fi
  done
fi
if [[ -z "$CONDA_EXE" ]]; then
  echo "Conda was not found. The existing Spyder/EAS setup suggests it is installed,"
  echo "but this helper could not locate it. Send the diagnostic ZIP rather than improvising."
  exit 2
fi
CONDA_BASE="$($CONDA_EXE info --base)"
# shellcheck disable=SC1090
source "$CONDA_BASE/etc/profile.d/conda.sh"
echo "Conda: $CONDA_EXE"

# Use a dedicated environment so the scientific Spyder 'eas' environment is untouched.
# IMPORTANT macOS/conda detail: compiler activation hooks from cctools inspect variables
# such as AR that may be unset in the parent shell. Bash `set -u` would treat that as
# a fatal error. Temporarily disable nounset only while conda runs its activation hooks.
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo
  echo "[SETUP] Creating isolated conda environment '$ENV_NAME'. Internet is required."
  conda create -y -n "$ENV_NAME" -c conda-forge \
    python=3.11 cmake ninja 'conan>=2,<3' git coreutils sed compilers make pkg-config \
    pyarrow pyyaml pandas numpy scipy
else
  echo
  echo "[SETUP] Reusing existing isolated conda environment '$ENV_NAME'."
  echo "[SETUP] The previous successful package transaction is preserved; no reinstall is needed."
fi

echo "[SETUP] Activating '$ENV_NAME' with macOS-safe shell handling."
set +u
conda activate "$ENV_NAME"
set -u

# Fail early if the environment is incomplete.  R3 exposed one additional upstream
# CORSIKA build-time Python dependency: src/framework/core/make_input_file.py imports
# the Scikit-HEP `particle` package while generating particle-property headers.
# Pin to particle 0.25.1, the version documented for CORSIKA 8 in the current framework
# paper, rather than accepting a newer potentially incompatible major release.
for tool in python cmake conan git ninja make; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Required tool missing after activation: $tool"; exit 2; }
done
PARTICLE_OK=0
if python - <<'PYCHECK' >/dev/null 2>&1
from importlib.metadata import version
import particle
raise SystemExit(0 if version('particle') == '0.25.1' else 1)
PYCHECK
then
  PARTICLE_OK=1
fi
if (( PARTICLE_OK == 0 )); then
  echo "[SETUP] Installing pinned CORSIKA Python build dependency: particle==0.25.1"
  echo "[SETUP] Internet is required for this small one-time Python package install."
  python -m pip install --disable-pip-version-check --no-cache-dir "particle==0.25.1"
fi
python - <<'PYENV'
from importlib.metadata import version
import numpy, pandas, pyarrow, scipy, yaml, particle
assert version('particle') == '0.25.1', version('particle')
print('[SETUP] Python scientific/build dependencies import successfully.')
print('[SETUP] particle version:', version('particle'))
PYENV

# R14 minimal-paper build retains the proven compiler/disk/Python portability lineage. R2 exposed a known GNU m4 1.4.19 / GCC 15
# incompatibility: GCC 15 defaults C to gnu23, while the old gnulib bundled in m4
# misdetects that standard. The published workaround is to compile m4 under gnu17.
# Prefer the macOS system Apple Clang toolchain for C/C++; retain conda gfortran.
SHIM="$WORK/system_compiler_shims"
mkdir -p "$SHIM"
ln -sf /usr/bin/clang "$SHIM/clang"
ln -sf /usr/bin/clang "$SHIM/cc"
ln -sf /usr/bin/clang "$SHIM/gcc"
ln -sf /usr/bin/clang++ "$SHIM/clang++"
ln -sf /usr/bin/clang++ "$SHIM/c++"
ln -sf /usr/bin/clang++ "$SHIM/g++"
export PATH="$SHIM:$PATH"
export CC="$SHIM/clang"
export CXX="$SHIM/clang++"
export CFLAGS="${CFLAGS:-} -std=gnu17"
if command -v gfortran >/dev/null 2>&1; then export FC="$(command -v gfortran)"; fi
echo "[TOOLS] Forced C compiler: $($CC --version | head -1)"
echo "[TOOLS] Forced C++ compiler: $($CXX --version | head -1)"
echo "[TOOLS] C standard workaround: -std=gnu17"

# Clear only non-critical Conan working folders left by failed/finished source builds.
# Installed package artifacts remain in the cache and can still be reused.
echo "[DISK] Cleaning non-critical Conan source/build/download/temp cache."
conan cache clean "*" --source --build --download --temp >/dev/null 2>&1 || true
FREE_KB=$(df -Pk "$HOME/Desktop" | awk 'NR==2 {print $4}')
FREE_GB=$((FREE_KB/1024/1024))
echo "[DISK] Free space after safe Conan cache cleanup: approximately ${FREE_GB} GiB"
# If R3 already completed the expensive Conan dependency installation, its generated
# toolchain plus the persistent log provide a deterministic resume marker.  In that
# case the package artifacts remain installed and we do not need the original 30-GiB
# pre-install headroom again. Production-shower storage will be reconsidered only after
# the manuscript-driven Stage 02 audit.
DEPS_ALREADY_INSTALLED=0
if [[ -f "$SRC/conan_cmake/conan_toolchain.cmake" && -f "$LOG" ]] && grep -q "Install finished successfully" "$LOG"; then
  DEPS_ALREADY_INSTALLED=1
  echo "[RESUME] Previous Conan dependency installation is complete and reusable."
fi
if (( DEPS_ALREADY_INSTALLED == 0 && FREE_GB < 30 )); then
  echo "A fresh dependency build requires at least 30 GiB free after safe cache cleanup."
  echo "Current free space is approximately ${FREE_GB} GiB."
  echo "STOP here. Free more internal space or use an external drive; do not force the build."
  exit 2
fi
if (( DEPS_ALREADY_INSTALLED == 1 && FREE_GB < 18 )); then
  echo "The dependency stage is reusable, but less than 18 GiB remains for the CORSIKA compile/smoke stage."
  echo "STOP here rather than risking the internal SSD."
  exit 2
fi

# Regenerate Conan profiles under the system Apple Clang shim. This prevents the
# stale R2 clang-21/GCC-15 profile from controlling the new build.
conan profile detect --force >/dev/null
conan profile detect --name corsika8 --force >/dev/null
echo "[CONAN] Re-detected profiles under system Apple Clang."
conan profile show -pr corsika8 || true

# Be compatible with either form of the upstream macOS fix: some CORSIKA source
# revisions call GNU tools as readlink/sed, while proposed macOS patches call
# them greadlink/gsed. The conda packages provide GNU readlink and sed; create
# local aliases only when the g-prefixed names are absent.
if ! command -v greadlink >/dev/null 2>&1 && command -v readlink >/dev/null 2>&1; then
  ln -sf "$(command -v readlink)" "$CONDA_PREFIX/bin/greadlink"
fi
if ! command -v gsed >/dev/null 2>&1 && command -v sed >/dev/null 2>&1; then
  ln -sf "$(command -v sed)" "$CONDA_PREFIX/bin/gsed"
fi

mkdir -p "$WORK"

echo
printf '[TOOLS] Python: '; python --version
printf '[TOOLS] CMake: '; cmake --version | head -1
printf '[TOOLS] Conan: '; conan --version
printf '[TOOLS] Git: '; git --version
printf '[TOOLS] C++: '; "${CXX:-c++}" --version | head -1 || true
printf '[TOOLS] Fortran: '; "${FC:-gfortran}" --version | head -1 || true

if [[ ! -d "$SRC/.git" ]]; then
  echo
  echo "[SOURCE] Cloning official CORSIKA 8 recursively. Internet is required."
  git clone --recurse-submodules "$CORSIKA_URL" "$SRC"
else
  echo
  if git -C "$SRC" cat-file -e "$CORSIKA_COMMIT^{commit}" 2>/dev/null; then
    echo "[SOURCE] Existing checkout already contains the locked commit; no network fetch needed."
  else
    echo "[SOURCE] Locked commit absent locally; refreshing refs (internet required)."
    git -C "$SRC" fetch --all --tags --prune
  fi
fi

echo "[SOURCE] Checking out the locked commit."
git -C "$SRC" checkout --detach "$CORSIKA_COMMIT"
# R14 reproducibility hardening: discard tracked edits left by earlier failed
# portability attempts, but deliberately KEEP untracked conan_cmake artifacts
# so the completed dependency transaction remains reusable.
git -C "$SRC" reset --hard "$CORSIKA_COMMIT"
git -C "$SRC" submodule sync --recursive
git -C "$SRC" submodule update --init --recursive
ACTUAL_COMMIT="$(git -C "$SRC" rev-parse HEAD)"
if [[ "$ACTUAL_COMMIT" != "$CORSIKA_COMMIT" ]]; then
  echo "Commit mismatch: $ACTUAL_COMMIT"
  exit 3
fi


# ---------------------------------------------------------------------------
# R14 PAPER-MINIMAL SOURCE LAYER
# ---------------------------------------------------------------------------
SUPPORT_CPP="$ROOT/support/c8_pentagon_gamma_minimal.cpp"
ACTUAL_CPP_SHA="$(shasum -a 256 "$SUPPORT_CPP" | awk '{print $1}')"
if [[ "$ACTUAL_CPP_SHA" != "$CUSTOM_CPP_SHA" ]]; then
  echo "Bundled minimal paper application failed its SHA-256 check."
  exit 3
fi
cp "$SUPPORT_CPP" "$SRC/applications/c8_air_shower.cpp"
echo "[PAPER-MINIMAL] Installed gamma-only application source as applications/c8_air_shower.cpp"

# Verify from the ACTUAL locked checkout—not from external documentation—that
# the SIBYLL Decay object already defaults to handling all decays.  R13's only
# unresolved final-link symbol was the redundant setHandleAllDecay() call.
# We fail closed if the locked header differs from that assumption.
SIBYLL_DECAY_HEADER="$SRC/corsika/modules/sibyll/Decay.hpp"
python - "$SIBYLL_DECAY_HEADER" "$SRC/applications/c8_air_shower.cpp" "$WORK/R14_SIBYLL_DECAY_DEFAULT_AUDIT.json" <<'PYR14DECAY'
from pathlib import Path
import hashlib,json,re,sys
hdr=Path(sys.argv[1]); app=Path(sys.argv[2]); out=Path(sys.argv[3])
if not hdr.is_file(): raise SystemExit(f"missing locked SIBYLL decay header: {hdr}")
hs=hdr.read_text()
aps=app.read_text()
rx=re.compile(r'\bbool\s+handleAllDecays_\s*=\s*true\s*;')
matches=rx.findall(hs)
if len(matches) != 1:
    raise SystemExit(
        "R14 decay-default guard: expected exactly one locked-source "
        "'bool handleAllDecays_ = true;' declaration; found "+str(len(matches))
    )
if 'setHandleAllDecay(' in aps:
    raise SystemExit("R14 minimal app still calls redundant setHandleAllDecay(); refusing to build.")
out.write_text(json.dumps({
    'scope':'locked-source SIBYLL decay default audit',
    'scientific_settings_changed':False,
    'header':str(hdr),
    'header_sha256':hashlib.sha256(hs.encode()).hexdigest(),
    'verified_default':'handleAllDecays_ = true',
    'minimal_app_calls_setHandleAllDecay':False
},indent=2,sort_keys=True)+'\n')
print('[PAPER-MINIMAL] Locked SIBYLL decay default verified: handleAllDecays_=true; redundant setter omitted.')
PYR14DECAY

# Apple libc++ portability in the locked CORSIKA logging formatter.  This is the
# exact final-translation-unit blocker exposed by R10; it changes no logging content
# or physics, only the vendor-private type spelling.
SPDLOG_SPECIAL="$SRC/corsika/detail/framework/core/SpdlogSpecializations.inl"
python - "$SPDLOG_SPECIAL" "$WORK/R14_SPDLOG_PORTABILITY_PATCH.json" <<'PYR14SPD'
from pathlib import Path
import hashlib,json,sys
p=Path(sys.argv[1]); out=Path(sys.argv[2])
if not p.is_file(): raise SystemExit(f"missing {p}")
s=p.read_text(); old='std::_Put_time<char>'; new='decltype(std::put_time(nullptr, ""))'
if s.count(old) != 1:
    raise SystemExit(f'R14 logging guard: expected exactly one {old!r}; found {s.count(old)}')
before=hashlib.sha256(s.encode()).hexdigest(); s=s.replace(old,new,1); p.write_text(s)
out.write_text(json.dumps({'scope':'Apple libc++ std::put_time formatter portability only','scientific_settings_changed':False,'path':str(p),'sha256_before':before,'sha256_after':hashlib.sha256(s.encode()).hexdigest()},indent=2,sort_keys=True)+'\n')
print('[MACOS-SOURCE] Portable std::put_time formatter type applied.')
PYR14SPD

# Preserve the proven macOS Conan libc++ correction.
python - "$SRC/conan-install.sh" <<'PYR14CONAN'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text(); s2=s.replace('compiler.libcxx=libstdc++11','compiler.libcxx=libc++')
if s2 != s: p.write_text(s2); print('[SOURCE] Applied macOS libc++ patch to conan-install.sh')
else: print('[SOURCE] No libstdc++11 profile line required patching at this commit.')
PYR14CONAN
chmod +x "$SRC/conan-install.sh" || true
[[ -f "$SRC/corsika-cmake.sh" ]] && chmod +x "$SRC/corsika-cmake.sh" || true

# The R12/R14 minimal-paper architecture: the standard CORSIKA8 umbrella target links every available
# interaction package into c8_air_shower. For this manuscript baseline we detach
# ONLY the explicitly unused module families from CORSIKA8. Their source remains
# untouched and available for a later, separately justified cross-check package.
#
# Retained baseline physics:
#   SIBYLL 2.3d (high-energy hadrons + decay), UrQMD (low-energy hadrons),
#   SOPHIA (resonance photonuclear), PROPOSAL (EM transport), core/output stack.
# Detached from this baseline executable:
#   Pythia8, TAUOLA, QGSJet-II, QGSJet-III, EPOS-LHC, EPOS-LHC-R, CONEX.
python - "$SRC" "$WORK/R14_MODULE_PRUNE.json" <<'PYR14PRUNE'
from pathlib import Path
import hashlib,json,re,sys
src=Path(sys.argv[1]); out=Path(sys.argv[2])
banned=('pythia','tauola','qgsjetii','qgsjetiii','epos-lhc','epos-lhcr','conex')
# Match complete CMake attachment commands, including multiline formatting.
rx=re.compile(r'(?ms)^(?P<indent>[ \t]*)(?P<cmd>add_dependencies|target_link_libraries)\s*\(\s*CORSIKA8\b.*?\)\s*$')
records=[]
for q in sorted(src.rglob('CMakeLists.txt')):
    try: before=q.read_text()
    except (OSError,UnicodeDecodeError): continue
    rel=str(q.relative_to(src)).lower(); changed=0
    pieces=[]; pos=0
    for m in rx.finditer(before):
        block=m.group(0); low=block.lower().replace('-','')
        belongs=any(f'modules/{b}' in rel for b in banned)
        names_banned=any(b.replace('-','') in low for b in banned)
        if not (belongs or names_banned): continue
        pieces.append(before[pos:m.start()])
        pieces.append('\n'.join('# R14 PAPER-MINIMAL DETACHED: '+line for line in block.splitlines()))
        pos=m.end(); changed+=1
    if changed:
        pieces.append(before[pos:]); after=''.join(pieces); q.write_text(after)
        records.append({'path':str(q),'commands_detached':changed,'sha256_before':hashlib.sha256(before.encode()).hexdigest(),'sha256_after':hashlib.sha256(after.encode()).hexdigest()})
remaining=[]
for q in sorted(src.rglob('CMakeLists.txt')):
    try: text=q.read_text()
    except (OSError,UnicodeDecodeError): continue
    rel=str(q.relative_to(src)).lower()
    for m in rx.finditer(text):
        block=m.group(0); low=block.lower().replace('-','')
        if any(f'modules/{b}' in rel for b in banned) or any(b.replace('-','') in low for b in banned):
            remaining.append(f'{q}: {block[:180]!r}')
if remaining:
    raise SystemExit('R14 module-prune verification failed; active banned attachment remains: '+repr(remaining[:20]))
if not records:
    raise SystemExit('R14 module-prune guard found no CORSIKA8 attachments to detach; refusing to guess.')
out.write_text(json.dumps({'scope':'paper-minimal CORSIKA8 dependency pruning','scientific_baseline_retained':['SIBYLL-2.3d','UrQMD','SOPHIA','PROPOSAL'],'detached_families':list(banned),'records':records},indent=2,sort_keys=True)+'\n')
print(f'[PAPER-MINIMAL] Detached unused CORSIKA8 module attachments in {len(records)} CMake file(s).')
for r in records: print('   ',r['path'], 'commands=',r['commands_detached'])
PYR14PRUNE

# SIBYLL is retained and its framework needs the known Darwin guard around GNU ld's
# ELF-only --exclude-libs option. Apply it only to SIBYLL; unused modules are not built.
python - "$SRC/modules/sibyll/CMakeLists.txt" "$WORK/R14_SIBYLL_LINK_PATCH.json" <<'PYR14SIB'
from pathlib import Path
import hashlib,json,re,sys
p=Path(sys.argv[1]); out=Path(sys.argv[2]); s=p.read_text()
marker='endif() # R14 Darwin: GNU/ELF-only --exclude-libs'
if marker in s:
    n=1; s2=s
else:
    rx=re.compile(r'^(?P<i>[ \t]*)(?P<l>target_link_options\(SIBYLL23d PRIVATE ["\']LINKER:--exclude-libs,ALL["\']\)[^\n]*)$',re.M)
    def repl(m):
        i=m.group('i'); l=m.group('l')
        return f'{i}if(NOT APPLE)\n{i}  {l.lstrip()}\n{i}{marker}'
    s2,n=rx.subn(repl,s)
    if n != 1: raise SystemExit(f'R14 SIBYLL linker guard expected one unguarded directive; found {n}')
    p.write_text(s2)
out.write_text(json.dumps({'scope':'SIBYLL framework Darwin linker portability','scientific_settings_changed':False,'path':str(p),'guarded_directives':n,'sha256':hashlib.sha256(s2.encode()).hexdigest()},indent=2,sort_keys=True)+'\n')
print('[MACOS-SOURCE] SIBYLL GNU/ELF-only --exclude-libs guarded on Apple.')
PYR14SIB

mkdir -p "$BUILD"
if (( DEPS_ALREADY_INSTALLED == 1 )); then
  echo "[RESUME] Skipping the already-successful Conan dependency installation."
  echo "[RESUME] Reusing installed package artifacts and generated Conan toolchain."
else
  echo "[BUILD] Installing pinned CORSIKA dependencies through the source tree's Conan recipe."
  ( cd "$BUILD"; "$SRC/conan-install.sh" --source-directory "$SRC" --release-with-debug )
fi

echo "[DISK] Cleaning non-critical Conan working folders before minimal CORSIKA compile."
conan cache clean "*" --source --build --download --temp >/dev/null 2>&1 || true
FREE_KB=$(df -Pk "$HOME/Desktop" | awk 'NR==2 {print $4}')
FREE_GB=$((FREE_KB/1024/1024))
echo "[DISK] Free before minimal CORSIKA compile: approximately ${FREE_GB} GiB"
if (( FREE_GB < 15 )); then
  echo "Less than 15 GiB remains; stopping safely."
  exit 2
fi

echo
echo "[BUILD] Configuring CORSIKA for the paper-minimal target."
( cd "$BUILD"; if [[ -x "$SRC/corsika-cmake.sh" ]]; then "$SRC/corsika-cmake.sh"; else cmake -S "$SRC" -B "$BUILD" -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo; fi )

# Verify the generated final application link is genuinely minimal BEFORE compiling.
LINKTXT="$BUILD/applications/CMakeFiles/c8_air_shower.dir/link.txt"
FLAGSTXT="$BUILD/applications/CMakeFiles/c8_air_shower.dir/flags.make"
[[ -f "$LINKTXT" && -f "$FLAGSTXT" ]] || { echo "Generated c8_air_shower link/flags files missing."; exit 3; }
python - "$LINKTXT" <<'PYR14LINKAUDIT'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text().lower()
banned=['pythia','tauola','qgsjetii','qgsjetiii','epos-lhc','epos_lhc','epos-lhcr','epos_lhcr','conex']
hits=[x for x in banned if x in s]
if hits: raise SystemExit(f'R14 minimal-link audit failed: unused module names remain in final link command: {hits}')
required=['sibyll23d','sophia','urqmd','proposal']
missing=[x for x in required if x not in s]
if missing: raise SystemExit(f'R14 minimal-link audit failed: required baseline component(s) absent from link command: {missing}')
if '--exclude-libs' in s or '-xlinker all' in s:
    raise SystemExit('R14 minimal-link audit failed: GNU/orphan linker token remains.')
print('[PAPER-MINIMAL] Generated link audit PASSED: SIBYLL + SOPHIA + UrQMD + PROPOSAL retained; broad model families absent.')
PYR14LINKAUDIT

# Fast logging/header preflight before invoking the dependency build.
PROBE_DIR="$WORK/macos_compile_probes"; mkdir -p "$PROBE_DIR"
cat > "$PROBE_DIR/minimal_logging_probe.cpp" <<'CPP'
#include <corsika/framework/core/Logging.hpp>
#include <iomanip>
#include <ctime>
int main(){ std::tm t{}; CORSIKA_LOG_INFO("{}", std::put_time(&t, "%Y")); return 0; }
CPP
python - "$FLAGSTXT" "$CXX" "$PROBE_DIR/minimal_logging_probe.cpp" "$PROBE_DIR/minimal_logging_probe.out" <<'PYR14PROBE'
from pathlib import Path
import shlex,subprocess,sys
flags=Path(sys.argv[1]).read_text().splitlines(); cxx=sys.argv[2]; src=sys.argv[3]; out=sys.argv[4]
vals={}
for line in flags:
    if ' = ' in line:
        k,v=line.split(' = ',1); vals[k.strip()]=v.strip()
cmd=[cxx]+shlex.split(vals.get('CXX_DEFINES',''))+shlex.split(vals.get('CXX_INCLUDES',''))+shlex.split(vals.get('CXX_FLAGS',''))+['-fsyntax-only',src]
r=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
Path(out).write_text('$ '+' '.join(shlex.quote(x) for x in cmd)+'\n'+r.stdout)
if r.returncode: raise SystemExit('R14 logging preflight failed; see '+out)
print('[MACOS-PREFLIGHT] Core logging/std::put_time syntax probe passed.')
PYR14PROBE

# R12 showed that a small header probe is not enough: application-level template/API
# errors can otherwise appear only when the target reaches 100%.  Build only the tiny
# generated-header/code-table prerequisites, then compile the COMPLETE minimal application
# with the exact generated CMake C++ flags in -fsyntax-only mode.  This changes no objects
# and runs before the expensive retained libraries are rebuilt.
JOBS=2
APP_PREFLIGHT="$WORK/application_preflight"; mkdir -p "$APP_PREFLIGHT"
echo "[APP-PREFLIGHT] Building only lightweight generated-header prerequisites."
for tgt in GenParticlesHeaders GenMediaProperties SourceDirLinkSib SourceDirLinkSoph; do
  cmake --build "$BUILD" --target "$tgt" --parallel "$JOBS" >/dev/null
done
python - "$FLAGSTXT" "$CXX" "$SRC/applications/c8_air_shower.cpp" "$APP_PREFLIGHT/c8_air_shower_full_syntax.out" <<'PYR14FULLAPP'
from pathlib import Path
import shlex,subprocess,sys
flags=Path(sys.argv[1]).read_text().splitlines(); cxx=sys.argv[2]; src=sys.argv[3]; out=Path(sys.argv[4])
vals={}
for line in flags:
    if ' = ' in line:
        k,v=line.split(' = ',1); vals[k.strip()]=v.strip()
cmd=[cxx]+shlex.split(vals.get('CXX_DEFINES',''))+shlex.split(vals.get('CXX_INCLUDES',''))+shlex.split(vals.get('CXX_FLAGS',''))+['-fsyntax-only',src]
r=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
out.write_text('$ '+' '.join(shlex.quote(x) for x in cmd)+'\n'+r.stdout)
if r.returncode:
    tail='\n'.join(r.stdout.splitlines()[-80:])
    print(tail)
    raise SystemExit('R14 FULL application syntax preflight failed; no long CORSIKA target build was started. See '+str(out))
print('[APP-PREFLIGHT] FULL minimal c8_air_shower.cpp syntax/template check PASSED.')
PYR14FULLAPP

# The minimal dependency graph and complete application source are now validated.
# Build only the application target; CMake will build only retained transitive dependencies.
echo
echo "[BUILD] Building paper-minimal c8_air_shower with $JOBS workers."
cmake --build "$BUILD" --target c8_air_shower --parallel "$JOBS"

if [[ ! -x "$BINARY" ]]; then
  echo "Expected minimal executable missing: $BINARY"
  exit 4
fi


BINARY_SHA="$(shasum -a 256 "$BINARY" | awk '{print $1}')"
echo "[BUILD] Compile/link complete: $BINARY"
echo "[BUILD] Binary SHA-256 before runtime repair: $BINARY_SHA"
exit 0
