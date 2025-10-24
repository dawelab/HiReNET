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
# HiReNET: sharedHOR_more.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Identify and visualize shared higher-order repeat (HOR)
#   patterns across multiple chromosomes. This version always
#   runs the “_all.R” multi-chromosome workflow, providing
#   comprehensive cross-chromosomal HOR network analysis.
#
# Workflow:
#   1) S6_HOR_newpatt_all.R
#   2) S7_HOR_shared_pattern.R
#   3) S8_Shared_HOR_plot_cus_all.R
#
# Input:
#   --datadir <dir>       Directory containing *.blat.sub files
#                         (typically from compareConsensus outputs)
#   --outdir <dir>        Output directory (created if missing)
#   --chr <list>          Comma-separated chromosome list
#                         (e.g., "chr1,chr2,chr8")
#
# Optional:
#   --letter <dir>        Directory containing per-bin name↔letter
#                         mapping files (default: <outdir>/mergebin_string_outputs)
#
# Dependencies:
#   - Rscript
#   - awk, grep, sed, find, cp
#
# Output (under --outdir):
#   - shared_HOR_patterns_all/
#   - shared_HOR_summary_all.txt
#   - plots_all/ (multi-chromosome visualizations)
#
# Example:
#   HiReNET sharedHOR_more \
#       --datadir AthCEN178_compare_consensusHOR_all/blat_sub \
#       --outdir AthCEN178_shared_out_all \
#       --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs \
#       --chr "chr1,chr2,chr3,chr4,chr5"
#
# Notes:
#   - Always executes the *_all.R versions of shared HOR scripts.
#   - Filters the letters directory to include only specified chromosomes.
#   - Designed for whole-genome HOR comparisons across multiple arrays.
# ==============================================================

die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

usage() {
  cat <<'EOF'
============================================================
                   HiReNET: sharedHOR_more
============================================================

Description:
  Perform multi-chromosome shared HOR analysis and plotting.
  Runs the “_all.R” workflow to identify shared HOR patterns
  across multiple chromosomes simultaneously.

------------------------------------------------------------
Usage:
  HiReNET sharedHOR_more \
      --datadir <path/to/blat_sub_dir> \
      --outdir <path/to/output_dir> \
      --chr <csv_list> [--letter <path/to/letters_dir>]

Required arguments:
  --datadir <dir>       Directory containing *.blat.sub files
  --outdir <dir>        Output directory (created if missing)
  --chr <csv>           Comma-separated chromosome list
                        (e.g., "chr1,chr2,chr3")

Optional arguments:
  --letter <dir>        Directory with per-bin name↔letter files  
                        (default: <outdir>/mergebin_string_outputs)
  -h, --help            Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • shared_HOR_patterns_all/  
  • shared_HOR_summary_all.txt  
  • plots_all/ (multi-chromosome HOR visualizations)

Example:
  HiReNET sharedHOR_more \
      --datadir AthCEN178_compare_consensusHOR_all/blat_sub \
      --outdir AthCEN178_shared_out_all \
      --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs \
      --chr "chr1,chr2,chr3,chr4,chr5"

============================================================
EOF
}

# Backward compatibility (for legacy usage)
print_help() { usage; }

# ---------- Parse args ----------
[[ $# -eq 0 ]] && { usage; exit 1; }

DATADIR=""
OUTDIR=""
LETTERDIR=""
CHR_FILTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datadir|--datadir=*)
      if [[ "$1" == *=* ]]; then DATADIR="${1#*=}"; else DATADIR="${2:-}"; shift; fi
      shift;;
    --outdir|--outdir=*)
      if [[ "$1" == *=* ]]; then OUTDIR="${1#*=}"; else OUTDIR="${2:-}"; shift; fi
      shift;;
    --letter|--letter=*)
      if [[ "$1" == *=* ]]; then LETTERDIR="${1#*=}"; else LETTERDIR="${2:-}"; shift; fi
      shift;;
    --chr|--chr=*)
      if [[ "$1" == *=* ]]; then CHR_FILTER="${1#*=}"; else CHR_FILTER="${2:-}"; shift; fi
      shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

[[ -n "$DATADIR"    ]] || die "Missing --datadir"
[[ -n "$OUTDIR"     ]] || die "Missing --outdir"
[[ -n "$CHR_FILTER" ]] || die "Missing --chr (comma-separated list required)"

# ---------- Deps ----------
need Rscript
need awk; need grep; need sed; need find; need cp

# ---------- Paths ----------
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

# ---------- Sanity ----------
[[ -d "$DATADIR"   ]] || die "--datadir not found: $DATADIR"
num_sub=$(find "$DATADIR" -maxdepth 1 -type f -name '*.blat.sub' | wc -l || true)
[[ "$num_sub" -gt 0 ]] || die "No *.blat.sub files in: $DATADIR"
[[ -d "$LETTERDIR" ]] || die "--letter directory not found: $LETTERDIR"

# ---------- Filter letters by --chr ----------
TMP_LETTERDIR="$OUTDIR/filtered_letters"
mkdir -p "$TMP_LETTERDIR"
IFS=',' read -ra CHRS <<< "$CHR_FILTER"
echo "[*] Filtering letters for chromosomes: ${CHRS[*]}"
found_any=0
for chr in "${CHRS[@]}"; do
  shopt -s nullglob
  matches=("$LETTERDIR"/${chr}*)
  if (( ${#matches[@]} == 0 )); then
    echo "WARN: No files matching ${chr}* in $LETTERDIR"
  else
    cp "${matches[@]}" "$TMP_LETTERDIR"/
    found_any=1
  fi
done
if (( found_any == 1 )); then
  LETTERDIR="$TMP_LETTERDIR"
else
  echo "WARN: No matching letter files copied; continuing with original LETTERDIR."
fi

echo "Params:"
echo "  datadir   : $DATADIR"
echo "  outdir    : $OUTDIR"
echo "  letters   : $LETTERDIR"
echo "  chr list  : ${CHRS[*]}"
echo

# ---------- Run multi-chr workflow (always *_all.R) ----------
echo "[*] S6_HOR_newpatt_all.R ..."
Rscript R/S6_HOR_newpatt_all.R \
  --input "$DATADIR" \
  --outdir "$OUTDIR"

echo "[*] S7_HOR_shared_pattern.R ..."
Rscript R/S7_HOR_shared_pattern.R \
  --outdir "$OUTDIR" \
  --letters "$LETTERDIR"

echo "[*] S8_Shared_HOR_plot_cus_all.R ..."
Rscript R/S8_Shared_HOR_plot_cus_all.R \
  --input "$DATADIR" \
  --outdir "$OUTDIR"

echo "[✓] Done. Outputs in: $OUTDIR"
