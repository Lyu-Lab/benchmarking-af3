#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda not found in PATH." >&2
  exit 1
fi

ARCH="$(uname -m)"
echo "[INFO] Detected architecture: $ARCH"

create_env() {
  local yml="$1"
  local name
  name="$(grep '^name:' "$yml" | awk '{print $2}')"
  if conda env list | awk '{print $1}' | grep -qx "$name"; then
    echo "[SKIP] $name already exists (remove with: conda env remove -n $name)"
    return 0
  fi
  conda env create -f "$yml"
}

echo "[INFO] Creating py3d (native)..."
create_env "$SCRIPT_DIR/envs/py3d_osx.yml"

if [[ "$ARCH" == "arm64" ]]; then
  echo "[INFO] Creating dock37_py27 under Rosetta (CONDA_SUBDIR=osx-64)..."
  CONDA_SUBDIR=osx-64 create_env "$SCRIPT_DIR/envs/dock37_py27_osx.yml"
else
  echo "[INFO] Creating dock37_py27 (native Intel Mac)..."
  create_env "$SCRIPT_DIR/envs/dock37_py27_osx.yml"
fi

cat <<EOF

[INFO] Environments ready. Before running the pipeline:

  export DOCKBASE=/path/to/DOCK-3.7-trunk

Then from enrichment/:

  conda activate py3d
  python split_running_sum.py --input-dir ./example_running_sum --out-dir ./example_running_sum
  bash compute_auc_from_splits.sh --root ./example_running_sum/docking_score

EOF
