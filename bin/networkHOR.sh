#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# networkHOR.sh
# Run HOR network (merged-bin) evaluation and summaries:
#   1) S3_network_summary_HOR_mergedbin.R
#   2) S4_HORpatterns.R
#   3) S5_Plot_HORmonomersize_freq.R
#
# Usage:
#   ./networkHOR.sh \
#     --blatsub DIR \
#     --bins bins_combined.txt \
#     --coor monomer_bed_inbin.txt \
#     --outdir OUTDIR \
#     [--default-thr 0.35]
#
# Notes:
#   - --blatsub should contain *.blat.sub files named like:
#       chrX_ARRAYSTART_ARRAYEND_BINSTART_BINEND.blat.sub
#     (or prefixed by all_)
#   - --bins must include columns: chr, start, end, array, threshold
#   - --coor is the combined BED-like table:
#       chr array_start array_end monomer_start monomer_end bin_start bin_end
# =========================================================

die(){ echo "ERROR: $*" >&2; exit 1; }

print_help() {
  sed -n '1,120p' "$0" | sed -n 's/^# \{0,1\}//p' | sed -n '1,120p'
}

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

# Load R (if you use modules on the cluster)
if command -v module >/dev/null 2>&1; then
  module load R/4.3.1-foss-2022a || true
fi
command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH."

echo "[*] blatsub : $BLATSUB"
echo "[*] bins    : $BINS"
echo "[*] coor    : $COOR"
echo "[*] outdir  : $OUTDIR"
echo "[*] default thr : $DEFAULT_THR"

# ---------- Step 3: Network re-eval over merged bins ----------
echo "[*] S3_network_summary_HOR_mergedbin.R ..."
Rscript classpredictionR/S3_network_summary_HOR_mergedbin.R \
  --input "$BLATSUB" \
  --bins "$BINS" \
  --outdir "$OUTDIR" \
  --default-thr "$DEFAULT_THR"


# The S3 script produces these two files in $OUTDIR:
S3_LABELS="$OUTDIR/label_HOR_re-eval_clusters.txt"
S3_REEVAL="$OUTDIR/HOR_re-eval_clusters.txt"

[[ -s "$S3_LABELS" ]] || die "Expected labels file not found: $S3_LABELS"
[[ -s "$S3_REEVAL" ]] || echo "[!] Warning: $S3_REEVAL not found (continuing)"

## ---------- Step 4: HOR patterns ----------
echo "[*] S4_HORpatterns.R ..."
Rscript classpredictionR/S4_HORpatterns.R \
  --groups-dir "$OUTDIR/out_csv" \
  --reval "$S3_REEVAL" \
  --coor "$COOR" \
  --outdir "$OUTDIR"

# ---------- Step 5: Monomer-size frequency plot ----------
echo "[*] S5_Plot_HORmonomersize_freq.R ..."
Rscript classpredictionR/S5_Plot_HORmonomersize_freq.R  --hor "$OUTDIR/HOR_coor_bed.txt" --outdir $OUTDIR/fre_plots

echo "[✓] Done. Outputs in: $OUTDIR"