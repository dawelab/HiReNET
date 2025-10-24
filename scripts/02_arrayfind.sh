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
# HiReNET: arrayfind.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Identify tandem-repeat arrays by aligning a consensus sequence
#   against a genome FASTA, merging nearby hits into array intervals,
#   and extracting their coordinates and sequences.
#
# Dependencies:
#   - seqkit
#   - blastn
#   - bedtools
#
# Input:
#   Genome FASTA, consensus FASTA
#
# Output:
#   Directory containing merged array coordinates and extracted sequences.
#
# Usage example:
#   HiReNET arrayfind -g Ey15.fasta \
#                     -c AthCEN178_consensus.fasta \
#                     -o AthCEN178_arrayout \
#                     -p AthCEN178
# ==============================================================

# ---- defaults ----
CHR_CSV=""
CHR_FILE=""
MIN_HIT_LEN=30
MERGE_DIST=10000
MIN_ARRAY_LEN=5000   # reserved for compatibility

# ---- helpers ----
die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

usage() {
  cat <<'EOF'
============================================================
                     HiReNET: arrayfind
============================================================

Description:
  Identify tandem repeat arrays across chromosomes by aligning
  a consensus repeat sequence to a genome FASTA. Arrays are 
  merged based on distance thresholds and filtered by length.

------------------------------------------------------------
Usage:
  HiReNET arrayfind \
      -g <genome.fasta> \
      -c <consensus.fasta> \
      -o <outdir> \
      -p <prefix> [options]

Required arguments:
  -g   Genome FASTA file
  -c   Consensus repeat FASTA file
  -o   Output directory (will be freshly created)
  -p   Prefix for output files

Optional arguments:
  --chr "chr1,chr2,chr3"     Comma-separated chromosome list to analyze
  --chr-file <file>          File containing chromosome names (one per line)
  --min-hit-len <int>        Minimum BLAT hit length to retain (default: 30)
  --merge-dist <int>         Merge nearby hits within this distance (default: 10000)
  --min-array-len <int>      Minimum array length (default: 5000, reserved)
  -h, --help                 Show this help message and exit

------------------------------------------------------------
Notes:
  • Automatically detects chr*.fa / chr*.fasta if no --chr or --chr-file is provided.
  • Deletes any existing output directory before starting (fresh run).
  • Ensures consistent BLAT filtering and merging across all chromosomes.

Example:
  HiReNET arrayfind \
      -g Ey15.fasta \
      -c AthCEN178_consensus.fasta \
      -o AthCEN178_arrayout \
      -p AthCEN178

============================================================
EOF
}

# Backward compatibility (if older scripts call print_help)
print_help() { usage; }

# ---- parse args ----
if [[ $# -eq 0 ]]; then print_help; exit 0; fi

GENOME=""
CONSENSUS=""
OUTBASE=""
PREFIX=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--genome)       GENOME="${2:-}"; shift 2;;
    -c|--consensus)    CONSENSUS="${2:-}"; shift 2;;
    -o|--outdir)       OUTBASE="${2:-}"; shift 2;;
    -p|--prefix)       PREFIX="${2:-}"; shift 2;;
    --chr)             CHR_CSV="${2:-}"; shift 2;;
    --chr-file)        CHR_FILE="${2:-}"; shift 2;;
    --min-hit-len)     MIN_HIT_LEN="${2:-}"; shift 2;;
    --merge-dist)      MERGE_DIST="${2:-}"; shift 2;;
    --min-array-len)   MIN_ARRAY_LEN="${2:-}"; shift 2;;
    -h|--help)         print_help; exit 0;;
    *) die "Unknown argument: $1 (use -h)";;
  esac
done

[[ -z "$GENOME"    ]] && die "Provide -g/--genome"
[[ -z "$CONSENSUS" ]] && die "Provide -c/--consensus"
[[ -z "$OUTBASE"   ]] && die "Provide -o/--outdir"
[[ -z "$PREFIX"    ]] && die "Provide -p/--prefix"

# ---- convert to absolute paths (Linux) ----
GENOME="$(readlink -f "$GENOME")"
CONSENSUS="$(readlink -f "$CONSENSUS")"
OUTDIR="$(readlink -f "${OUTBASE}")"

# ---- clean previous output ----
if [[ -d "$OUTDIR" ]]; then
  echo "[!] Removing existing output directory: $OUTDIR"
  rm -rf "$OUTDIR"
fi
mkdir -p "$OUTDIR"

# ---- check dependencies ----
need seqkit; need blastn; need bedtools
need awk; need cut; need sort; need basename

# ---- output structure ----
SPLIT_DIR="$OUTDIR/split_seq"
BLAST_DIR="$OUTDIR/blast"
ARRAY_DIR="$OUTDIR"
mkdir -p "$SPLIT_DIR" "$BLAST_DIR" "$ARRAY_DIR"

# ---- split genome ----
echo "[*] Splitting genome with seqkit..."
seqkit split -i --by-id-prefix "" -O "$SPLIT_DIR" "$GENOME"
if [[ -d "$SPLIT_DIR/split" ]]; then
  SPLIT_DIR="$SPLIT_DIR/split"
fi

# ---- determine chromosome files ----
declare -a CHR_PATHS
if [[ -n "$CHR_FILE" ]]; then
  mapfile -t wanted < <(grep -v '^[[:space:]]*$' "$CHR_FILE")
  for id in "${wanted[@]}"; do
    if   [[ -f "$SPLIT_DIR/${id}.fasta" ]]; then CHR_PATHS+=("$SPLIT_DIR/${id}.fasta")
    elif [[ -f "$SPLIT_DIR/${id}.fa"    ]]; then CHR_PATHS+=("$SPLIT_DIR/${id}.fa")
    else echo "WARN: requested $id not found in $SPLIT_DIR (skip)"; fi
  done
elif [[ -n "$CHR_CSV" ]]; then
  IFS=',' read -r -a wanted <<< "$CHR_CSV"
  for id in "${wanted[@]}"; do
    if   [[ -f "$SPLIT_DIR/${id}.fasta" ]]; then CHR_PATHS+=("$SPLIT_DIR/${id}.fasta")
    elif [[ -f "$SPLIT_DIR/${id}.fa"    ]]; then CHR_PATHS+=("$SPLIT_DIR/${id}.fa")
    else echo "WARN: requested $id not found in $SPLIT_DIR (skip)"; fi
  done
else
  shopt -s nullglob
  files=( "$SPLIT_DIR"/chr*.fa "$SPLIT_DIR"/chr*.fasta )
  shopt -u nullglob
  if ((${#files[@]} == 0)); then
    die "No chr*.fa or chr*.fasta files found in $SPLIT_DIR"
  fi

  # Sort by natural chromosome order
  mapfile -t CHR_PATHS < <(
    printf '%s\n' "${files[@]}" \
    | awk '{
        path=$0
        key=path
        sub(/^.*\//,"",key)
        sub(/\.(fa|fasta)$/,"",key)
        print key "\t" path
      }' \
    | sort -V -k1,1 \
    | cut -f2
  )
fi

((${#CHR_PATHS[@]})) || die "No chromosome FASTAs found/chosen."
echo "[*] Chromosomes to process:"
printf '    %s\n' "${CHR_PATHS[@]##*/}"

# ---- run BLAST and generate arrays ----
for chr_path in "${CHR_PATHS[@]}"; do
  fname="$(basename "$chr_path")"
  chr_base="${fname%.fa}"
  chr_base="${chr_base%.fasta}"

  blast_txt="$BLAST_DIR/${chr_base}_${PREFIX}_blast.txt"
  merged_bed="$BLAST_DIR/${chr_base}_${PREFIX}_blast_filter_merge.bed"
  array_bed="$ARRAY_DIR/${chr_base}_${PREFIX}_array.bed"
  array_fa="$ARRAY_DIR/${chr_base}_${PREFIX}_array.fa"

  echo "[*] BLAST: $chr_base vs $CONSENSUS"
  blastn -subject "$chr_path" -query "$CONSENSUS" -outfmt 6 -max_target_seqs 5000000 > "$blast_txt"

  echo "[*] Filter+merge BLAST hits → BED ($chr_base)"
  awk -v MINLEN="$MIN_HIT_LEN" '$4>=MINLEN' "$blast_txt" \
    | cut -f2,9,10 \
    | awk '{if($3<$2){start=$3-1; end=$2}else{start=$2-1; end=$3}; print $1"\t"start"\t"end}' \
    | bedtools sort -i - \
    | bedtools merge -d "$MERGE_DIST" -i - -c 1 -o count \
    | awk '{if($4>=10) print $0"\t"$3-$2}' \
    | awk '{if($5/($3-$2)>.1) print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4}' > "$merged_bed"

  cp "$merged_bed" "$array_bed"

  if [[ -s "$merged_bed" ]]; then
    echo "[*] Extract array FASTA ($chr_base)"
    bedtools getfasta -fi "$chr_path" -bed "$merged_bed" -fo "$array_fa"
  else
    echo "[!] No merged intervals for $chr_base (skipping getfasta)."
    : > "$array_fa"
  fi
done

echo "[✓] Done."
echo "Output base:      $OUTDIR"
echo "  Split FASTAs:   $SPLIT_DIR"
echo "  BLAST results:  $BLAST_DIR/*_${PREFIX}_blast.txt"
echo "  Array BEDs:     $ARRAY_DIR/*_${PREFIX}_array.bed"
echo "  Array FASTAs:   $ARRAY_DIR/*_${PREFIX}_array.fa"
