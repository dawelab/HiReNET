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
# HiReNET: sharedHOR.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Identify and visualize shared higher-order repeat (HOR)
#   patterns across one or more chromosomes. Integrates output
#   from previous HOR consensus and classification steps to
#   detect shared sequence motifs and relationships.
#
# Workflow:
#   1) Run S1_HOR_newpatt.R
#   2) Run S2_HOR_shared_pattern.R
#   3) Run S3_Shared_HOR_plot_cus.R
#
# Input:
#   --datadir <dir>       Directory containing *.blat.sub files
#                         (typically from compareConsensus)
#   --outdir <dir>        Output directory (created if missing)
#
# Optional:
#   --letter <dir>        Directory containing per-bin name↔letter
#                         mapping files (default: <outdir>/mergebin_string_outputs)
#   --chr <list>          Comma-separated list of chromosomes
#                         (e.g., "chr1" or "chr1,chr2")
#
# Dependencies:
#   - Rscript
#   - awk, grep, sed, find, cp
#
# Output (under --outdir):
#   - shared_HOR_patterns/
#   - shared_HOR_summary.txt
#   - plots/ (custom R visualizations)
#
# Example:
#   HiReNET sharedHOR \
#       --datadir AthCEN178_compare_consensusHOR_chr1/blat_sub \
#       --outdir AthCEN178_shared_out_chr1 \
#       --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs \
#       --chr "chr1"
#
# Notes:
#   - If --chr is provided, only matching letter files (e.g., chr1*)
#     are copied and analyzed.
#   - R scripts are expected under R/sharedhorR/.
# ==============================================================

die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing dependency in PATH: $1"; }

usage() {
  cat <<'EOF'
============================================================
                     HiReNET: sharedHOR
============================================================

Description:
  Detect and visualize shared HOR patterns across one or more
  chromosomes. Uses outputs from consensus and network steps
  to identify cross-chromosomal HOR relationships.

------------------------------------------------------------
Usage:
  HiReNET sharedHOR \
      --datadir <path/to/blat_sub_dir> \
      --outdir <path/to/output_dir> [options]

Required arguments:
  --datadir <dir>       Directory containing *.blat.sub files
  --outdir <dir>        Output directory (created if missing)

Optional arguments:
  --letter <dir>        Directory containing per-bin name↔letter files  
                        (default: <outdir>/mergebin_string_outputs)
  --chr <list>          Comma-separated chromosome list  
                        (e.g., "chr1" or "chr1,chr2")
  --plotv <V>  Plot versions, V1, V2, V3
  -h, --help            Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • shared_HOR_patterns/  
  • shared_HOR_summary.txt  
  • plots/ (visualizations from R scripts)

Example:
  HiReNET sharedHOR \
      --datadir AthCEN178_compare_consensusHOR_chr1/blat_sub \
      --outdir AthCEN178_shared_out_chr1 \
      --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs \
      --chr "chr1"

============================================================
EOF
}

# Backward compatibility
print_help() { usage; }
# ---------- Parse args ----------
[[ $# -eq 0 ]] && { usage; exit 1; }

DATADIR=""
OUTDIR=""
LETTERDIR=""
CHR_FILTER=""
PLOTV="" 

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
    --plotv)        PLOTV="${2:-}"; shift 2;;
    --plotv=*)      PLOTV="${1#*=}"; shift;;
    -h|--help)   usage; exit 0;;
    *)           die "Unknown argument: $1";;
  esac
done

[[ -n "$DATADIR" ]] || die "Missing --datadir"
[[ -n "$OUTDIR"  ]] || die "Missing --outdir"

# ---------- Normalize paths & deps ----------
need Rscript
need find
need awk
need grep
need sed
need cp

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

# ---------- Sanity checks ----------
[[ -d "$DATADIR"   ]] || die "--datadir not found: $DATADIR"
num_sub=$(find "$DATADIR" -maxdepth 1 -type f -name '*.blat.sub' | wc -l || true)
[[ "$num_sub" -gt 0 ]] || die "No *.blat.sub files in: $DATADIR"
[[ -d "$LETTERDIR" ]] || die "--letter directory not found: $LETTERDIR"

# ---------- Apply optional --chr filter to letters ----------
ACTIVE_LETTERDIR="$LETTERDIR"
if [[ -n "$CHR_FILTER" ]]; then
  TMP_LETTERDIR="$OUTDIR/filtered_letters"
  mkdir -p "$TMP_LETTERDIR"
  IFS=',' read -ra CHRS <<< "$CHR_FILTER"
  echo "[*] Filtering letters for chromosomes: ${CHRS[*]}"
  found_any=0
  shopt -s nullglob
  for chr in "${CHRS[@]}"; do
    matches=( "$LETTERDIR"/${chr}* )
    if (( ${#matches[@]} == 0 )); then
      echo "WARN: No files matching ${chr}* in $LETTERDIR"
    else
      cp "${matches[@]}" "$TMP_LETTERDIR"/
      found_any=1
    fi
  done
  if (( found_any == 1 )); then
    ACTIVE_LETTERDIR="$TMP_LETTERDIR"
  else
    echo "WARN: No matching letter files were copied; using original LETTERDIR."
  fi
fi

echo "Params:"
echo "  datadir   : $DATADIR"
echo "  outdir    : $OUTDIR"
echo "  letters   : $ACTIVE_LETTERDIR"
echo "  chr filter: ${CHR_FILTER:-<none>}"
echo

# ---------- Run ----------
echo "[*] S6_HOR_newpatt.R ..."
Rscript ${HIRENET_ROOT}/R/S6_HOR_newpatt.R \
  --input "$DATADIR" \
  --outdir "$OUTDIR"

echo "[*] S7_HOR_shared_pattern.R ..."
Rscript ${HIRENET_ROOT}/R/S7_HOR_shared_pattern.R \
  --outdir "$OUTDIR" \
  --letters "$ACTIVE_LETTERDIR"

echo "[*] S8_Shared_HOR_plot_cus2.R ..."
Rscript "${HIRENET_ROOT}/R/S8_Shared_HOR_plot_cus2.R" \
  --outdir "$OUTDIR" \
  --plotv "$PLOTV"

echo "[✓] Done. Outputs in: $OUTDIR"
