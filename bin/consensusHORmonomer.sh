#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# consensusHORmonomer.sh
# - Build per-class consensus (clustalo + EMBOSS cons)
# - Concatenate per-chromosome consensus
# - Split per-threshold (robust for 0.90 vs 0.9)
# - BLAT: within-chr and within-threshold
#
# Usage:
#   ./consensusHORmonomer.sh \
#     --outdir network_mergebin_consensus \
#     --chroms "chr1,chr2" \
#     --thresholds "0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99" \
#     --threads 4 \
#     --minseq 2
#
# Expects OUTDIR/group_HOR_monomers_output from arrangeHORmonomer.sh
# Requires: clustalo, EMBOSS 'cons', bioawk, blat, awk, grep
# ------------------------------------------------------------

OUTDIR=""
CHROMS="chr1,chr2,chr3,chr4,chr5"
THRESHOLDS="0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99"
MINSEQ=2
THREADS=2
SKIP_MODULES=0

usage() {
  cat <<'USAGE'
consensusHORmonomer.sh
Required:
  --outdir         Output base used by arrangeHORmonomer.sh (contains group_HOR_monomers_output)

Optional:
  --chroms         Comma list (default: "chr1,chr2,chr3,chr4,chr5")
  --thresholds     Comma list (default: "0.91...,0.99")
  --minseq         Min sequences to build a consensus (default: 2)
  --threads        Threads for clustalo (default: 2)
  --skip-modules   Skip `ml` loads
  -h|--help        Show this help

Example:
  ./consensusHORmonomer.sh --outdir network_mergebin_consensus --chroms "chr1,chr2" --threads 4
USAGE
  exit 1
}

[[ $# -eq 0 ]] && usage
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)        OUTDIR="${2:-}"; shift 2;;
    --outdir=*)      OUTDIR="${1#*=}"; shift;;
    --chroms)        CHROMS="${2:-}"; shift 2;;
    --chroms=*)      CHROMS="${1#*=}"; shift;;
    --thresholds)    THRESHOLDS="${2:-}"; shift 2;;
    --thresholds=*)  THRESHOLDS="${1#*=}"; shift;;
    --minseq)        MINSEQ="${2:-}"; shift 2;;
    --minseq=*)      MINSEQ="${1#*=}"; shift;;
    --threads)       THREADS="${2:-}"; shift 2;;
    --threads=*)     THREADS="${1#*=}"; shift;;
    --skip-modules)  SKIP_MODULES=1; shift;;
    -h|--help)       usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -z "$OUTDIR" ]] && { echo "ERROR: --outdir is required"; usage; }

# load modules (optional)
if [[ $SKIP_MODULES -eq 0 ]] && command -v ml >/dev/null 2>&1; then
  ml Clustal-Omega/1.2.4-GCC-12.3.0 || true
  ml EMBOSS/6.6.0-foss-2023a || true
  ml BLAT/3.7-GCC-12.3.0 || true
  ml bioawk/1.0-GCC-12.3.0 || true
fi

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing <$1>" >&2; exit 2; }; }
need awk; need grep; need clustalo; need cons; need bioawk; need blat

OUTDIR="$(readlink -f "$OUTDIR")"

# Inputs/outputs inside OUTDIR
out_lists="$OUTDIR/group_HOR_monomers_output"
[[ -d "$out_lists" ]] || { echo "ERROR: missing $out_lists (run arrangeHORmonomer.sh first)"; exit 3; }

cons_dir="$OUTDIR/group_HOR_monomers_consensus"
all_cons_dir="$OUTDIR/all_recluster_consensus_monomer"
op_group_dir="$all_cons_dir/op_group"
blat_chr_out="$all_cons_dir/blat_HOR_con_output"
blat_chr_sub="$all_cons_dir/blat_HOR_con_output_sub"
blat_thr_out="$all_cons_dir/blat_HOR_con_output2"
blat_thr_sub="$all_cons_dir/blat_HOR_con_output_sub2"

mkdir -p "$cons_dir" "$all_cons_dir" "$op_group_dir" \
         "$blat_chr_out" "$blat_chr_sub" "$blat_thr_out" "$blat_thr_sub"

IFS=',' read -r -a CHR_ARR <<< "$CHROMS"
IFS=',' read -r -a THR_ARR <<< "$THRESHOLDS"

echo "Params:"
echo "  OUTDIR      = $OUTDIR"
echo "  CHROMS      = $CHROMS"
echo "  THRESHOLDS  = $THRESHOLDS"
echo "  MINSEQ      = $MINSEQ"
echo "  THREADS     = $THREADS"
echo

# --------- 2) Build consensus per class (letter) ----------
shopt -s nullglob
for d in "$out_lists"/*; do
  [[ -d "$d" ]] || continue
  dirbase="$(basename "$d")"           # e.g., chr1_..._10000
  outbin="$cons_dir/$dirbase"
  mkdir -p "$outbin"

  for file in "$d"/*.fa; do
    [[ -e "$file" ]] || continue
    base="$(basename "$file" .fa)"     # e.g., chr1_..._10000_B_0.92
    aln="$outbin/${base}_clust.fasta"
    consout="$outbin/${base}.clust_CONS.fa"

    nseq=$(grep -c '^>' "$file" || true)
    if [[ "${nseq:-0}" -lt "$MINSEQ" ]]; then
      echo "WARN: <$file> has <$nseq> sequences; skipping consensus." >&2
      continue
    fi

    [[ -s "$aln" ]] || clustalo -i "$file" -o "$aln" --force --threads "$THREADS"

    # Header = exact basename (no prefixing)
    label="$base"
    [[ -s "$consout" ]] || cons -sequence "$aln" -name "$label" -outseq "$consout"
  done
done

# --------- 3) Concatenate consensus per chromosome ----------
for chr in "${CHR_ARR[@]}"; do
  outfa="$all_cons_dir/${chr}_clust.cons.fa"
  : > "$outfa"
  for dd in "$cons_dir"/${chr}*; do
    [[ -d "$dd" ]] || continue
    cat "$dd"/*.clust_CONS.fa 2>/dev/null || true
  done >> "$outfa"
  echo "Wrote: $outfa"
done

# --------- 4) Split per-threshold (robust to 0.90 vs 0.9) ----------
for chr in "${CHR_ARR[@]}"; do
  src="$all_cons_dir/${chr}_clust.cons.fa"
  [[ -s "$src" ]] || { echo "WARN: missing $src"; continue; }
  for th in "${THR_ARR[@]}"; do
    th_out="${th/0./}"        # 0.90 -> 90
    thA="$th"
    thB="${th%0}"             # 0.90 -> 0.9; 0.97 -> 0.97 (unchanged)

    dst="$op_group_dir/${chr}_clust.cons.${th_out}.fa"
    bioawk -c fastx -v thA="$thA" -v thB="$thB" '
      index($name, thA) || index($name, thB) { print ">"$name"\n"$seq }
    ' "$src" > "$dst"
  done
done

for th in "${THR_ARR[@]}"; do
  th_out="${th/0./}"
  out="$op_group_dir/clust.cons.${th_out}.fa"
  cat "$op_group_dir"/chr*_clust.cons."$th_out".fa > "$out" 2>/dev/null || true
  echo "[*] Wrote $out"
done



echo "[✓] consensusHORmonomer: done"
echo "  Per-bin consensus         : $cons_dir/"
echo "  Per-chr consensus FASTAs  : $all_cons_dir/*_clust.cons.fa"
echo "  Split per-threshold       : $op_group_dir/"
