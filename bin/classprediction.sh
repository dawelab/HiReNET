#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# classprediction.sh
# Run the HOR classification pipeline:
#   1) S1_network_unmerge_bin_{noplot,plot}.R  (chosen via --plot)
#   2) S2_smoothing_bin.R
#
# Usage:
#   ./classprediction.sh --blatsub DIR --outdir DIR --prefix NAME [--bin N] [--plot]
#   or:
#   ./classprediction.sh --blat DIR   --outdir DIR --prefix NAME [--bin N] [--plot]
#
# Notes:
#   - If --blat points to .blat files, they are converted to .blat.sub first.
#   - --plot selects the plotting S1 script; omit to run the no-plot S1.
#   - --bin controls S2 only (default: 10000).
# =========================================================

die() { echo "ERROR: $*" >&2; exit 1; }

print_help() {
  cat <<'EOF'
Usage:
  classprediction.sh (--blat DIR | --blatsub DIR) --outdir DIR --prefix NAME [--bin N] [--plot]

Options:
  --blat DIR        Directory with *.blat files (will be converted to *.blat.sub)
  --blatsub DIR     Directory with *.blat.sub files (preferred)
  --outdir DIR      Output directory (will be created)
  --prefix NAME     Prefix used by S1/S2 outputs (e.g., AthCEN178)
  --bin N           Bin size for S2_smoothing_bin.R (default: 10000)
  --plot            Use plotting S1 (S1_network_unmerge_bin_plot.R). Without this, uses no-plot S1.
  -h, --help        Show this help
EOF
}

# -------- args --------
[[ $# -eq 0 ]] && { print_help; exit 1; }

BLATDIR=""
OUTDIR=""
PREFIX=""
MODE=""      # "blat" or "blatsub"
PLOT=0       # 0=no, 1=yes
BIN=10000    # default bin size for S2

while [[ $# -gt 0 ]]; do
  case "$1" in
    --blat)       BLATDIR="${2:-}"; MODE="blat"; shift 2;;
    --blatsub)    BLATDIR="${2:-}"; MODE="blatsub"; shift 2;;
    --outdir)     OUTDIR="${2:-}"; shift 2;;
    --prefix)     PREFIX="${2:-}"; shift 2;;
    --bin)        BIN="${2:-}"; shift 2;;
    --plot)       PLOT=1; shift 1;;
    -h|--help)    print_help; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "${BLATDIR}" ]] && die "Provide --blat DIR or --blatsub DIR"
[[ -z "${OUTDIR}"  ]] && die "Provide --outdir DIR"
[[ -z "${PREFIX}"  ]] && die "Provide --prefix NAME"

BLATDIR="$(readlink -f "$BLATDIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"

# -------- R availability / module --------
if command -v module >/dev/null 2>&1; then
  module load R/4.3.1-foss-2022a || true
fi
command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH."

echo "[*] BLATDIR : $BLATDIR"
echo "[*] OUTDIR  : $OUTDIR"
echo "[*] PREFIX  : $PREFIX"
echo "[*] MODE    : ${MODE:-blatsub}"
echo "[*] PLOT    : $([[ $PLOT -eq 1 ]] && echo yes || echo no)"
echo "[*] BIN     : $BIN"

# -------- if needed, convert .blat -> .blat.sub --------
if [[ "$MODE" == "blat" ]]; then
  echo "[*] Converting .blat files to .blat.sub..."
  SUBDIR="$OUTDIR/blat_sub"
  mkdir -p "$SUBDIR"
  shopt -s nullglob
  found=0
  for blat in "$BLATDIR"/*.blat; do
    found=1
    name=$(basename "$blat" .blat)
    awk '$1 ~ /^[0-9]+$/' "$blat" \
      | awk 'NR>1' \
      | awk '{print $1"\t"$10"\t"$11"\t"$14"\t"$15}' \
      > "$SUBDIR/${name}.blat.sub"
  done
  [[ $found -eq 1 ]] || die "No .blat files found in $BLATDIR"
  BLATDIR="$SUBDIR"
fi

# -------- choose S1 script based on --plot --------
if [[ $PLOT -eq 1 ]]; then
  S1_SCRIPT="classpredictionR/S1_network_unmerge_bin_plot.R"
else
  S1_SCRIPT="classpredictionR/S1_network_unmerge_bin_noplot.R"
fi
[[ -f "$S1_SCRIPT" ]] || die "Missing S1 script: $S1_SCRIPT"

# -------- S1 --------
echo "[*] Running $(basename "$S1_SCRIPT") ..."
Rscript "$S1_SCRIPT" \
        --input "$BLATDIR" \
        --outdir "$OUTDIR" \
        --prefix "$PREFIX" \
        --lda-model classpredictionR/lda_model.rds

# -------- S2 --------
echo "[*] Running S2_smoothing_bin.R ..."
Rscript classpredictionR/S2_smoothing_bin.R \
        --input "$OUTDIR/${PREFIX}_bin_class.txt" \
        --outdir "$OUTDIR" \
        --prefix "$PREFIX" \
        --bin "$BIN"

echo "[✓] Done. Outputs in $OUTDIR"