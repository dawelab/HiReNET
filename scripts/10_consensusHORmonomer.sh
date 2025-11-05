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
# HiReNET: consensusHORmonomer.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Build per-class, per-chromosome consensus sequences for HOR 
#   monomer groups. Each consensus is derived from the outputs of 
#   arrangeHORmonomer.sh using multiple sequence alignment 
#   (Clustal Omega) followed by consensus generation (EMBOSS cons).
#
# Workflow:
#   1) Generate consensus FASTAs for each HOR group using Clustal Omega + cons
#   2) Concatenate per-chromosome consensus sequences
#   3) Split consensus sets by identity threshold (e.g., 0.90–0.99)
#
# Input:
#   --outdir <dir>        Directory containing group_HOR_monomers_output/
#                         (produced by arrangeHORmonomer.sh)
#
# Optional parameters:
#   --chroms <list>       Comma-separated chromosome list
#                         (default: chr1–chr10)
#   --thresholds <list>   Comma-separated identity thresholds
#                         (default: 0.90–0.99)
#   --minseq <int>        Minimum number of sequences to build consensus
#                         (default: 2)
#   --threads <int>       Threads for Clustal Omega alignment (default: 2)
#
# Dependencies:
#   - Clustal Omega (clustalo)
#   - EMBOSS cons
#   - bioawk, awk, grep, sed
#
# Output (under --outdir):
#   - consensus_per_chr/*.fasta
#   - all_recluster_consensus_monomer/
#       ├── chrX_threshold0.90_consensus.fa
#       ├── chrX_threshold0.91_consensus.fa
#       ├── ...
#
# Example:
#   HiReNET consensusHORmonomer \
#       --outdir AthCEN178_network_mergebin_consensus \
#       --chroms "chr1,chr2" \
#       --thresholds "0.90,0.95,0.99" \
#       --threads 8
#
# Notes:
#   - Requires output from HiReNET arrangeHORmonomer (step 09).
#   - Automatically handles threshold-based splitting for HOR robustness.
# ==============================================================

OUTDIR=""
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10"
THRESHOLDS="0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99"
MINSEQ=2
THREADS=2

usage() {
  cat <<'EOF'
============================================================
                HiReNET: consensusHORmonomer
============================================================

Description:
  Build per-chromosome HOR monomer consensus sequences 
  using Clustal Omega and EMBOSS cons.  
  Generates consensus sets across multiple identity thresholds.

------------------------------------------------------------
Usage:
  HiReNET consensusHORmonomer \
      --outdir <path/to/output_dir> [options]

Required arguments:
  --outdir <dir>          Directory containing group_HOR_monomers_output/
                          (produced by arrangeHORmonomer)

Optional arguments:
  --chroms <list>         Comma-separated list of chromosomes
                          (default: "chr1,chr2,...,chr10")
  --thresholds <list>     Comma-separated identity thresholds
                          (default: "0.90,0.91,...,0.99")
  --minseq <int>          Minimum number of sequences to build consensus (default: 2)
  --threads <int>         Threads for Clustal Omega (default: 2)
  -h, --help              Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • consensus_per_chr/*.fasta  
  • all_recluster_consensus_monomer/
      ├── chrX_threshold0.90_consensus.fa  
      ├── chrX_threshold0.95_consensus.fa  
      ├── chrX_threshold0.99_consensus.fa  

Example:
  HiReNET consensusHORmonomer \
      --outdir AthCEN178_network_mergebin_consensus \
      --chroms "chr1,chr2" \
      --threads 4

============================================================
EOF
}

# Backward compatibility (for older usage)
print_help() { usage; }

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
    -h|--help)       usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -z "$OUTDIR" ]] && { echo "ERROR: --outdir is required"; usage; }

# --- deps (PATH-based) ---
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing <$1>" >&2; exit 2; }; }
need awk; need grep; need sed
need clustalo
need cons
need bioawk

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

# --------- 2) Build consensus per class (letter, original) ----------
#shopt -s nullglob
#for d in "$out_lists"/*; do
#  [[ -d "$d" ]] || continue
#  dirbase="$(basename "$d")"           # e.g., chr1_..._10000
#  outbin="$cons_dir/$dirbase"
#  mkdir -p "$outbin"
#
#  for file in "$d"/*.fa; do
#    [[ -e "$file" ]] || continue
#    base="$(basename "$file" .fa)"     # e.g., chr1_..._10000_B_0.92
#    aln="$outbin/${base}_clust.fasta"
#    consout="$outbin/${base}.clust_CONS.fa"
#
#    nseq=$(grep -c '^>' "$file" || true)
#    if [[ "${nseq:-0}" -lt "$MINSEQ" ]]; then
#      echo "WARN: <$file> has <$nseq> sequences; skipping consensus." >&2
#      continue
#    fi
#
#    [[ -s "$aln" ]] || clustalo -i "$file" -o "$aln" --force --threads "$THREADS"
#
#    # Header = exact basename (no prefixing)
#    label="$base"
#    [[ -s "$consout" ]] || cons -sequence "$aln" -name "$label" -outseq "$consout"
#  done
#done

# --------- 2) Build consensus per class (letter, change 1) ----------
#shopt -s nullglob
#for d in "$out_lists"/*; do
#  [[ -d "$d" ]] || continue
#  dirbase="$(basename "$d")"            # e.g., chr1_..._10000
#  chr_prefix="${dirbase%%_*}"           # "chr1" from "chr1_..."
#
#  # --- keep only requested chromosomes ---
#  keep=0
#  for want in "${CHR_ARR[@]}"; do
#    want_strip="${want//[[:space:]]/}"  # be robust to stray spaces
#    if [[ "$chr_prefix" == "$want_strip" ]]; then
#      keep=1
#      break
#    fi
#  done
#  [[ $keep -eq 1 ]] || continue
#
#  outbin="$cons_dir/$dirbase"
#  mkdir -p "$outbin"
#
#  for file in "$d"/*.fa; do
#    [[ -e "$file" ]] || continue
#    base="$(basename "$file" .fa)"     # e.g., chr1_..._10000_B_0.92
#    aln="$outbin/${base}_clust.fasta"
#    consout="$outbin/${base}.clust_CONS.fa"
#
#    nseq=$(grep -c '^>' "$file" || true)
#    if [[ "${nseq:-0}" -lt "$MINSEQ" ]]; then
#      echo "WARN: <$file> has <$nseq> sequences; skipping consensus." >&2
#      continue
#    fi
#
#    [[ -s "$aln"     ]] || clustalo -i "$file" -o "$aln" --force --threads "$THREADS"
#    [[ -s "$consout" ]] || cons -sequence "$aln" -name "$base" -outseq "$consout"
#  done
#done

# --------- 2) Build consensus per class (letter, change2) ----------
shopt -s nullglob

# Select directories only for specified chromosomes
chr_dirs=()
for chr in "${CHR_ARR[@]}"; do
  # Trim spaces and find all dirs beginning with chr_
  chr_clean="${chr//[[:space:]]/}"
  for d in "$OUTDIR/group_HOR_monomers_output/${chr_clean}_"*; do
    [[ -d "$d" ]] && chr_dirs+=("$d")
  done
done

if [[ ${#chr_dirs[@]} -eq 0 ]]; then
  echo "ERROR: No directories found matching requested chromosomes under $OUTDIR/group_HOR_monomers_output"
  exit 3
fi

for d in "${chr_dirs[@]}"; do
  dirbase="$(basename "$d")"            # e.g., chr1_14389697_14926924_14599767_14619767
  chr_prefix="${dirbase%%_*}"           # "chr1"
  echo "[*] Processing: $dirbase ($chr_prefix)"
  outbin="$cons_dir/$dirbase"
  mkdir -p "$outbin"

  for file in "$d"/*.fa; do
    [[ -e "$file" ]] || continue
    base="$(basename "$file" .fa)"
    aln="$outbin/${base}_clust.fasta"
    consout="$outbin/${base}.clust_CONS.fa"

    nseq=$(grep -c '^>' "$file" || true)
    if [[ "${nseq:-0}" -lt "$MINSEQ" ]]; then
      echo "WARN: <$file> has <$nseq> sequences; skipping consensus." >&2
      continue
    fi

    [[ -s "$aln"     ]] || clustalo -i "$file" -o "$aln" --force --threads "$THREADS"
    [[ -s "$consout" ]] || cons -sequence "$aln" -name "$base" -outseq "$consout"
  done
done
# --------- 3) Concatenate consensus per chromosome ----------
for chr in "${CHR_ARR[@]}"; do
  outfa="$all_cons_dir/${chr}_clust.cons.fa"
  : > "$outfa"
  for dd in "$cons_dir"/"${chr}_"*; do   # only dirs starting with 'chr_'
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
    th_out="${th/0./}"           # 0.90 -> 90 (filename suffix)
    thA="$th"                    # e.g. 0.90
    thB="${th%0}"                # e.g. 0.9  (handles files with one decimal)
    dst="$op_group_dir/${chr}_clust.cons.${th_out}.fa"

    bioawk -c fastx -v thA="$thA" -v thB="$thB" '
      {
        n = split($name, a, "_");
        thr = a[n];
        if (thr == thA || thr == thB) {
          print ">"$name"\n"$seq;
        }
      }
    ' "$src" > "$dst"
  done
done

# Merge all chromosomes per threshold
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
