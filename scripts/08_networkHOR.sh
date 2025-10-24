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
# HiReNET: networkHOR.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Perform HOR network evaluation, pattern detection, and 
#   frequency summarization across merged genomic bins. This stage 
#   integrates monomer similarity data to identify HOR structures 
#   and their monomer composition.
#
# Workflow:
#   1) S3_network_summary_HOR_mergedbin.R
#        - Re-evaluates network connectivity and HOR clusters
#   2) S4_HORpatterns.R
#        - Extracts monomer composition patterns for each HOR
#   3) S5_Plot_HORmonomersize_freq.R
#        - Summarizes monomer length and HOR frequency distributions
#
# Input:
#   --blatsub DIR      Directory containing *.blat.sub files
#                      (named as chrX_ARRAYSTART_ARRAYEND_BINSTART_BINEND.blat.sub)
#   --bins FILE        Combined bin summary table with columns:
#                      chr, start, end, array, threshold
#   --coor FILE        Monomer coordinate table:
#                      chr array_start array_end monomer_start monomer_end bin_start bin_end
#   --outdir DIR       Output directory (created if missing)
#
# Dependencies:
#   - Rscript (with required R packages)
#   - awk, sort
#
# Output (under --outdir):
#   - label_HOR_re-eval_clusters.txt
#   - HOR_re-eval_clusters.txt
#   - HOR_coor_bed.txt
#   - fre_plots/ (monomer frequency plots)
#
# Example:
#   HiReNET networkHOR \
#       --blatsub AthCEN178_comparemonomers/blat_output_sub \
#       --bins AthCEN178_classpred_out/AthCEN178_bin_class_smooth.txt \
#       --coor AthCEN178_arrangemonomer_10kb/AthCEN178_monomer_bed_inbin.txt \
#       --outdir AthCEN178_networkHOR_out \
#       --default-thr 0.35
#
# Notes:
#   - The default threshold (--default-thr) defines connectivity strength 
#     for HOR cluster detection (default: 0.35).
#   - Output files serve as input for HOR consensus construction (step 09).
# ==============================================================
die(){ echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat <<'EOF'
============================================================
                     HiReNET: networkHOR
============================================================

Description:
  Perform network-based HOR analysis across merged bins.  
  Identifies, summarizes, and visualizes higher-order repeat 
  clusters and their monomer composition.

------------------------------------------------------------
Usage:
  HiReNET networkHOR \
      --blatsub <path/to/blat_sub_dir> \
      --bins <path/to/bins_table.txt> \
      --coor <path/to/monomer_bed_inbin.txt> \
      --outdir <path/to/output_dir> [options]

Required arguments:
  --blatsub <dir>       Directory with *.blat.sub files
  --bins <file>         Combined bins summary (chr, start, end, array, threshold)
  --coor <file>         Monomer coordinate table
  --outdir <dir>        Output directory (created if missing)

Optional arguments:
  --default-thr <float> Default connectivity threshold (default: 0.35)
  -h, --help            Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • label_HOR_re-eval_clusters.txt  
  • HOR_re-eval_clusters.txt  
  • HOR_coor_bed.txt  
  • fre_plots/ (monomer frequency visualizations)

Example:
  HiReNET networkHOR \
      --blatsub AthCEN178_comparemonomers/blat_output_sub \
      --bins AthCEN178_classpred_out/AthCEN178_bin_class_smooth.txt \
      --coor AthCEN178_arrangemonomer_10kb/AthCEN178_monomer_bed_inbin.txt \
      --outdir AthCEN178_networkHOR_out \
      --default-thr 0.35

============================================================
EOF
}

# Backward compatibility (if older scripts still call print_help)
print_help() { usage; }

[[ $# -eq 0 ]] && { print_help; exit 1; }

BLATSUB=""
BINS=""
COOR=""
OUTDIR=""
DEFAULT_THR="0.35"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --blatsub)     BLATSUB="${2:-}"; shift 2;;
    --bins)        BINS="${2:-}"; shift 2;;
    --coor)        COOR="${2:-}"; shift 2;;
    --outdir)      OUTDIR="${2:-}"; shift 2;;
    --default-thr) DEFAULT_THR="${2:-}"; shift 2;;
    -h|--help)     print_help; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "$BLATSUB" ]] && die "Provide --blatsub"
[[ -z "$BINS"    ]] && die "Provide --bins"
[[ -z "$COOR"    ]] && die "Provide --coor"
[[ -z "$OUTDIR"  ]] && die "Provide --outdir"

BLATSUB="$(readlink -f "$BLATSUB")"
BINS="$(readlink -f "$BINS")"
COOR="$(readlink -f "$COOR")"
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"

# PATH-based check for R
command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH. Activate your env."

echo "[*] blatsub : $BLATSUB"
echo "[*] bins    : $BINS"
echo "[*] coor    : $COOR"
echo "[*] outdir  : $OUTDIR"
echo "[*] default thr : $DEFAULT_THR"

# ---------- Step 3: Network re-eval over merged bins ----------
echo "[*] S3_network_summary_HOR_mergedbin.R ..."
Rscript R/S3_network_summary_HOR_mergedbin.R \
  --input "$BLATSUB" \
  --bins "$BINS" \
  --outdir "$OUTDIR" \
  --default-thr "$DEFAULT_THR"

# Expected S3 outputs:
S3_LABELS="$OUTDIR/label_HOR_re-eval_clusters.txt"
S3_REEVAL="$OUTDIR/HOR_re-eval_clusters.txt"

[[ -s "$S3_LABELS" ]] || die "Expected labels file not found: $S3_LABELS"
[[ -s "$S3_REEVAL" ]] || echo "[!] Warning: $S3_REEVAL not found (continuing)"

# ---------- Step 4: HOR patterns ----------
echo "[*] S4_HORpatterns.R ..."
Rscript R/S4_HORpatterns.R \
  --groups-dir "$OUTDIR/out_csv" \
  --reval "$S3_REEVAL" \
  --coor "$COOR" \
  --outdir "$OUTDIR"

# ---------- Step 5: Monomer-size frequency plot ----------
echo "[*] S5_Plot_HORmonomersize_freq.R ..."
Rscript R/S5_Plot_HORmonomersize_freq.R \
  --hor "$OUTDIR/HOR_coor_bed.txt" \
  --outdir "$OUTDIR/fre_plots"

echo "[✓] Done. Outputs in: $OUTDIR"
