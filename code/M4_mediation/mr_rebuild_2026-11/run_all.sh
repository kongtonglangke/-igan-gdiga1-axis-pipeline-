#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Reproduce the rebuilt methylation -> IgAN MR layer (Figure 4 A-F, Table S5b-S5d).
#
#   ./run_all.sh                       # Rscript must be on PATH
#   RSCRIPT=/f/R-4.3.3/bin/x64/Rscript.exe ./run_all.sh
#
# Every path below is exported, never hard-coded inside the scripts; 00_config.R
# reads these variables and falls back to <MR_ROOT>/... when they are unset.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"

# R on Windows decodes `source()` paths with the local code page: a non-ASCII
# directory name makes every script fail with "cannot open the connection".
if printf '%s' "$HERE" | LC_ALL=C grep -q '[^ -~]'; then
  echo "ERROR: this path contains non-ASCII characters:" >&2
  echo "  $HERE" >&2
  echo "Clone the repository to an ASCII-only path (e.g. C:/repos/igan-pipeline) and re-run." >&2
  exit 1
fi

export MR_ROOT="$REPO/data/M4_mediation/mr_rebuild_2026-11"
export MR_DATA_DIR="$MR_ROOT/inputs"
export MQTL_DIR="$REPO/data/M4_mediation"      # the five cpg_full_assoc_<CpG>.tsv
export MR_FIG_DIR="$REPO/results/M4_mediation"
export EQTL_CACHE="${EQTL_CACHE:-$REPO/data/eQTL_cache}"

mkdir -p "$MR_ROOT/out" "$MR_FIG_DIR"

if [ ! -d "$MQTL_DIR" ]; then
  echo "MQTL_DIR not found: $MQTL_DIR" >&2
  exit 1
fi

for s in 01_instruments 02_mr_main 03_coloc 04_directionality \
         05_positive_control 06_reverse 07_simulation 08_figures; do
  printf '== %s ==\n' "$s"
  "$RSCRIPT" "$HERE/$s.R"
done

printf '\nDone. Tables: %s\n      Figures: %s\n' "$MR_ROOT/out" "$MR_FIG_DIR"
