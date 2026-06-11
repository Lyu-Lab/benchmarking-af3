#!/usr/bin/env bash
set -euo pipefail

# Ensure `conda activate` works in non-interactive shells by initializing
# conda for this script. This mirrors what `conda init` would do for an
# interactive shell but is safe to run in batch scripts.
if command -v conda >/dev/null 2>&1; then
  # Prefer the shell hook (conda >= 4.6). If that fails, fall back to
  # sourcing the profile.d script from the conda base prefix.
  eval "$(conda shell.bash hook 2>/dev/null)" || \
    source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" || true
fi

JOBS=1
DRYRUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="${2:-$PWD}"; shift 2 ;;
    --jobs|-j) JOBS="${2:-1}"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

# Resolve helper script paths relative to this script's location so the
# helpers are found even if the user runs this script from a different
# current working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONVERT_PY="$SCRIPT_DIR/convert.py"
RUNR_SH="$SCRIPT_DIR/runr.py"

if [[ -z "${DOCKBASE:-}" ]]; then
  echo "[ERROR] \$DOCKBASE is not set." >&2
  exit 2
fi

for f in "$CONVERT_PY" "$RUNR_SH"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing required file: $f" >&2; exit 2; }
done

is_highest_better_metric() {
  case "$1" in
    iptm_complex_af3|ptm_ligand_only_af3) return 0 ;; 
    *) return 1 ;;
  esac
}

process_one_dir() {
  local recp_dir="$1"
  local metric_dir="$(dirname "$recp_dir")"
  local metric="$(basename "$metric_dir")"
  local split_csv="$recp_dir/split.csv"
  local ligs="$recp_dir/ligands.name"
  local decs="$recp_dir/decoys.name"

  [[ -f "$split_csv" ]] || { echo "[SKIP] no split.csv: $recp_dir"; return 0; }
  [[ -s "$ligs" && -s "$decs" ]] || { echo "[SKIP] missing ligands/decoys: $recp_dir"; return 0; }

  echo "[INFO] >>> $metric / $(basename "$recp_dir")"

  if [[ "$DRYRUN" -eq 1 ]]; then
    if is_highest_better_metric "$metric"; then
      echo "  Would rank descending: (cd $recp_dir && conda activate dock37_py27 && conda activate py3d && python $RUNR_SH split.csv && conda deactivate)"
    else
      echo "  Would rank ascending:  (cd $recp_dir && conda activate dock37_py27 && conda activate py3d && python $CONVERT_PY split.csv && conda deactivate)"
    fi
    echo "  Would AUC/logAUC: (cd $recp_dir && python "$DOCKBASE/analysis/enrich.py" -i . -l ./ligands.name -d ./decoys.name)"
    echo "                     (cd $recp_dir && python "$DOCKBASE/analysis/plots.py"  -i . -l ./ligands.name -d ./decoys.name)"
    return 0
  fi

  if is_highest_better_metric "$metric"; then
    ( cd "$recp_dir" && conda activate dock37_py27 && conda activate py3d && python "$RUNR_SH" "split.csv" && conda deactivate )
  else
    ( cd "$recp_dir" && conda activate dock37_py27 && conda activate py3d && python "$CONVERT_PY" "split.csv" && conda deactivate )
  fi

  ( cd "$recp_dir" && conda activate dock37_py27 && python "$DOCKBASE/analysis/enrich.py" -i . -l ./ligands.name -d ./decoys.name )
  ( cd "$recp_dir" && conda activate dock37_py27 && python "$DOCKBASE/analysis/plots.py" -i . -l ./ligands.name -d ./decoys.name )
}

export CONVERT_PY RUNR_SH DRYRUN
export -f process_one_dir
export -f is_highest_better_metric


TARGET_DIRS=()
while IFS= read -r d; do
  TARGET_DIRS+=("$d")
done < <(find "$ROOT" -mindepth 2 -maxdepth 2 -type f -name "split.csv" -exec dirname {} \; | sort -u)
echo "[INFO] Found ${#TARGET_DIRS[@]} target directories under $ROOT"

if [[ "$DRYRUN" -eq 1 ]]; then
  for d in "${TARGET_DIRS[@]}"; do process_one_dir "$d"; done
  exit 0
fi

if [[ "$JOBS" -le 1 ]]; then
  for d in "${TARGET_DIRS[@]}"; do process_one_dir "$d"; done
else
  printf "%s\0" "${TARGET_DIRS[@]}" | xargs -0 -n1 -P "$JOBS" -I{} bash -c 'process_one_dir "$@"' _ {}
fi

echo "[DONE] All done."
