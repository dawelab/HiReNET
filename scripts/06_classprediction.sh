#!/usr/bin/env bash
set -euo pipefail

# --- find the HiReNET folder, but don't change working directory ---
if [[ -n "${HIRENET_ROOT:-}" && -d "${HIRENET_ROOT}" ]]; then
  ROOT="${HIRENET_ROOT}"
else
  ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
export HIRENET_ROOT="${ROOT}"
# ==============================================================
# HiReNET: classprediction.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Predict higher-order repeat (HOR) classes from BLAT-based 
#   monomer similarity networks. This step integrates BLAT output 
#   with R-based classification to identify HOR-rich regions and 
#   refine bin-level annotations.
#
# Workflow:
#   1) S1_network_unmerge_bin_{noplot,plot}.R
#        - Processes monomer similarity networks per bin
#        - Selects either plotting or non-plotting version via --plot
#   2) S2_smoothing_bin.R
#        - Smooths classification outputs across genomic bins
#
# Dependencies:
#   - Rscript (with required R libraries installed)
#   - awk, sed
#   - Optional: BLAT (if using --blat instead of --blatsub)
#
# Input:
#   --blat DIR        Directory with *.blat files (converted to .blat.sub)
#   --blatsub DIR     Directory with *.blat.sub files (preferred input)
#   --prefix NAME     Prefix for labeling all outputs
#
# Output (under --outdir):
#   - ${PREFIX}_bin_class.txt
#   - ${PREFIX}_bin_class_smooth.txt
#   - S1_* and S2_* R-generated visualizations (if --plot used)
#
# Example:
#   HiReNET classprediction \
#       --blatsub AthCEN178_comparemonomers/blat_output_sub \
#       --outdir  AthCEN178_classpred_out \
#       --prefix  AthCEN178 \
#       --bin 10000 \
#       --plot
#
# Notes:
#   - The --plot flag enables additional visualization in R (S1_plot mode).
#   - Both R scripts must be located in the R/ directory under HiReNET.
#   - Expected R outputs are used by rearrangemonomers.sh (next step).
# ==============================================================

die(){ echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat <<'EOF'
============================================================
                    HiReNET: classprediction
============================================================

Description:
  Predict higher-order repeat (HOR) classes based on BLAT-derived
  monomer similarity networks. Uses two R stages:
    1) S1_network_unmerge_bin_{noplot,plot}.R
    2) S2_smoothing_bin.R

------------------------------------------------------------
Usage:
  HiReNET classprediction \
      (--blat <dir> | --blatsub <dir>) \
      --outdir <dir> \
      --prefix <name> [options]

Required arguments:
  --blat <dir>       Directory with *.blat files (will be converted to *.blat.sub)
  --blatsub <dir>    Directory with *.blat.sub files (preferred)
  --outdir <dir>     Output directory (created if missing)
  --prefix <name>    Prefix for output files (e.g., AthCEN178)

Optional arguments:
  --bin <int>        Bin size for smoothing (default: 10000)
  --plot             Use plotting version of S1 (S1_network_unmerge_bin_plot.R)
  -h, --help         Show this help message and exit

------------------------------------------------------------
Notes:
  • If --blat is used, the script converts each .blat to .blat.sub automatically.  
  • The --plot option runs the plotting version of the R script; omit for no-plot mode.  
  • Ensure Rscript is in PATH and required R dependencies are installed.  

Example:
  HiReNET classprediction \
      --blatsub AthCEN178_comparemonomers/blat_output_sub \
      --outdir  AthCEN178_classpred_out \
      --prefix  AthCEN178 \
      --bin 10000 \
      --plot

============================================================
EOF
}

# Backward compatibility (if older scripts still call helpmsg or print_help)
print_help() { usage; }
helpmsg() { usage; }

# ---------- args ----------
[[ $# -eq 0 ]] && { usage; exit 1; }

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
    --plot)       PLOT=1; shift;;
    -h|--help)    usage; exit 0;;
    *)            die "Unknown argument: $1 (use -h for help)";;
  esac
done

# ---------- validation ----------
[[ -z "$BLATDIR" ]] && die "Provide --blat DIR or --blatsub DIR"
[[ -z "$OUTDIR"  ]] && die "Provide --outdir DIR"
[[ -z "$PREFIX"  ]] && die "Provide --prefix NAME"

# ---------- normalize paths ----------
BLATDIR="$(readlink -f "$BLATDIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"

command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH."

echo "[*] BLATDIR : $BLATDIR"
echo "[*] OUTDIR  : $OUTDIR"
echo "[*] PREFIX  : $PREFIX"
echo "[*] MODE    : ${MODE:-blatsub}"
echo "[*] PLOT    : $([[ $PLOT -eq 1 ]] && echo yes || echo no)"
echo "[*] BIN     : $BIN"

# convert .blat -> .blat.sub when needed
if [[ "$MODE" == "blat" ]]; then
  echo "[*] Converting .blat to .blat.sub ..."
  SUBDIR="$OUTDIR/blat_sub"
  mkdir -p "$SUBDIR"
  shopt -s nullglob
  cnt=0
  for blat in "$BLATDIR"/*.blat; do
    cnt=$((cnt+1))
    name=$(basename "$blat" .blat)
    awk '$1 ~ /^[0-9]+$/' "$blat" | awk 'NR>1' \
      | awk '{print $1"\t"$10"\t"$11"\t"$14"\t"$15}' \
      > "$SUBDIR/${name}.blat.sub"
  done
  [[ $cnt -gt 0 ]] || die "No .blat files found in $BLATDIR"
  BLATDIR="$SUBDIR"
fi

# ---- R scripts live directly under R/ now ----
S1_SCRIPT="R/S1_network_unmerge_bin_$([[ $PLOT -eq 1 ]] && echo plot || echo noplot).R"
S2_SCRIPT="R/S2_smoothing_bin.R"

# LDA model: prefer R/; fall back to scripts/ (matches your screenshots)
if   [[ -f "R/lda_model.rds" ]];        then LDA_MODEL="R/lda_model.rds"
elif [[ -f "scripts/lda_model.rds" ]];  then LDA_MODEL="scripts/lda_model.rds"
else die "Missing LDA model: put lda_model.rds in R/ (preferred) or scripts/."
fi

[[ -f "$S1_SCRIPT" ]] || die "Missing: $S1_SCRIPT"
[[ -f "$S2_SCRIPT" ]] || die "Missing: $S2_SCRIPT"

echo "[*] Running $(basename "$S1_SCRIPT") ..."
Rscript "$S1_SCRIPT" \
  --input  "$BLATDIR" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX" \
  --lda-model "$LDA_MODEL"

echo "[*] Running $(basename "$S2_SCRIPT") ..."
Rscript "$S2_SCRIPT" \
  --input  "$OUTDIR/${PREFIX}_bin_class.txt" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX" \
  --bin    "$BIN"

echo "[✓] Done. Outputs in $OUTDIR"
