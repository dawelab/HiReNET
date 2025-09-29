#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# sharedHOR.sh
# Run shared-HOR discovery/plotting with optional --chr filter
# =========================================================

die(){ echo "ERROR: $*" >&2; exit 1; }

usage(){
  cat <<'EOF'
Usage:
  sharedHOR.sh --datadir DIR --outdir DIR [--letter DIR] [--chr LIST]

Required:
  --datadir   Directory containing *.blat.sub files
  --outdir    Output directory (will be created)

Optional:
  --letter    Directory with per-bin nam↔letter files (default: <outdir>/mergebin_string_outputs)
  --chr       Comma-separated chromosome list.
              Example: --chr "chr1"   → only chr1* files from --letter
                       --chr "chr1,chr2" → chr1* and chr2* files

Notes:
  - Tries to load R via environment modules if available; otherwise uses PATH.
EOF
}

# ---------- Parse args ----------
[[ $# -eq 0 ]] && { usage; exit 1; }

DATADIR=""
OUTDIR=""
LETTERDIR=""
CHR_FILTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datadir)   DATADIR="${2:-}"; shift 2;;
    --datadir=*) DATADIR="${1#*=}"; shift;;
    --outdir)    OUTDIR="${2:-}"; shift 2;;
    --outdir=*)  OUTDIR="${1#*=}"; shift;;
    --letter)    LETTERDIR="${2:-}"; shift 2;;
    --letter=*)  LETTERDIR="${1#*=}"; shift;;
    --chr)       CHR_FILTER="${2:-}"; shift 2;;
    --chr=*)     CHR_FILTER="${1#*=}"; shift;;
    -h|--help)   usage; exit 0;;
    *)           die "Unknown argument: $1";;
  esac
done

[[ -n "$DATADIR" ]] || die "Missing --datadir"
[[ -n "$OUTDIR"  ]] || die "Missing --outdir"

# ---------- Normalize paths ----------
DATADIR="$(readlink -f "$DATADIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"

# default letters dir if not provided
if [[ -z "${LETTERDIR}" ]]; then
  LETTERDIR="$OUTDIR/mergebin_string_outputs"
fi
if [[ -d "$LETTERDIR" ]]; then
  LETTERDIR="$(readlink -f "$LETTERDIR")"
fi

# ---------- Apply chr filter ----------
TMP_LETTERDIR="$OUTDIR/filtered_letters"
if [[ -n "$CHR_FILTER" ]]; then
  mkdir -p "$TMP_LETTERDIR"
  IFS=',' read -ra CHRS <<< "$CHR_FILTER"
  echo "[*] Filtering letters for chromosomes: ${CHRS[*]}"
  for chr in "${CHRS[@]}"; do
    matches=$(ls "$LETTERDIR"/${chr}* 2>/dev/null || true)
    if [[ -z "$matches" ]]; then
      echo "WARN: No files matching ${chr}* in $LETTERDIR"
    else
      cp $matches "$TMP_LETTERDIR"/
    fi
  done
  LETTERDIR="$TMP_LETTERDIR"
fi

# ---------- Sanity checks ----------
[[ -d "$DATADIR" ]] || die "--datadir not found: $DATADIR"
num_sub=$(find "$DATADIR" -maxdepth 1 -type f -name '*.blat.sub' | wc -l || true)
[[ "$num_sub" -gt 0 ]] || die "No *.blat.sub files in: $DATADIR"
[[ -d "$LETTERDIR" ]] || die "--letter directory not found after filtering: $LETTERDIR"

# ---------- R availability ----------
if command -v module >/dev/null 2>&1; then
  module load R/4.3.1-foss-2022a || true
fi
command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH"

echo "Params:"
echo "  datadir   : $DATADIR"
echo "  outdir    : $OUTDIR"
echo "  letters   : $LETTERDIR"
echo "  chr filter: ${CHR_FILTER:-<none>}"
echo

# ---------- Run ----------
echo "[*] S1_HOR_newpatt.R ..."
Rscript sharedhorR/S1_HOR_newpatt.R \
  --input "$DATADIR" \
  --outdir "$OUTDIR"

echo "[*] S2_HOR_shared_pattern.R ..."
Rscript sharedhorR/S2_HOR_shared_pattern.R \
  --outdir "$OUTDIR" \
  --letters "$LETTERDIR"

echo "[*] S3_Shared_HOR_plot.R ..."
Rscript sharedhorR/S3_Shared_HOR_plot.R \
  --outdir "$OUTDIR"

echo "[✓] Done. Outputs in: $OUTDIR"