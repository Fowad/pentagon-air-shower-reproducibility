#!/bin/bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PENTAGON_CORSIKA_SOFTWARE_WORK="${PENTAGON_CORSIKA_SOFTWARE_WORK:-$HOME/Desktop/Pentagon_CORSIKA_Response_Stage1}"
# Prefer the already validated software build if present. Scientific outputs are NOT reused by the next launcher.
B="$PENTAGON_CORSIKA_SOFTWARE_WORK/corsika8-build/applications/c8_air_shower"
S="$PENTAGON_CORSIKA_SOFTWARE_WORK/corsika8-source/applications/c8_air_shower.cpp"
if [[ -x "$B" && -f "$S" ]]; then
  echo "[BOOTSTRAP] Existing CORSIKA software installation found; validating runtime instead of rebuilding software."
  bash "$ROOT/RUNTIME_AND_SMOKE.command"
  exit 0
fi
# On a genuinely clean Mac, install Miniforge if conda is absent, then perform the locked source build.
if ! command -v conda >/dev/null 2>&1 && [[ ! -x "$HOME/miniforge3/bin/conda" ]]; then
  echo "[BOOTSTRAP] Conda not found. Installing Miniforge3 x86_64 into ~/miniforge3 (internet required)."
  tmp="$(mktemp -t miniforge.XXXXXX.sh)"
  curl -L --fail --retry 3 -o "$tmp" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh"
  bash "$tmp" -b -p "$HOME/miniforge3"
  rm -f "$tmp"
fi
bash "$ROOT/BUILD_CORSIKA_FROM_CLEAN_SOURCE.command"
bash "$ROOT/RUNTIME_AND_SMOKE.command"
echo "[BOOTSTRAP] Locked CORSIKA software build/runtime is ready."
