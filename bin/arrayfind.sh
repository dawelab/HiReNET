#!/bin/bash
set -euo pipefail

# ======================================
# arrayfind.sh
# Find tandem-repeat arrays by blasting a consensus against a genome,
# merging nearby hits into array intervals, and extracting sequences.
#
# Example:
#   ./arrayfind.sh \
#     -g genome.fasta \
#     -c AthCEN178_consensus.fasta \
#     -o arrayout \
#     -p AthCEN178
#   # -> outputs in: AthCEN178_arrayout/
#
# Requires: seqkit, blastn, bedtools
# (Optional) modules: SeqKit, BLAST+, BEDTools
# ======================================

# ---- defaults ----
CHR_CSV=""
CHR_FILE=""
MIN_HIT_LEN=30
MERGE_DIST=100000
MIN_ARRAY_LEN=10000

# ---- helpers ----
die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

print_help() {
cat <<'EOF'
Usage:
  ./arrayfind.sh \
    -g /path/genome.fasta \
    -c /path/consensus.fasta \
    -o output_dir_base \
    -p PREFIX \
    [--chr "chr1,chr2,chr3"] \
    [--chr-file chr_list.txt] \
    [--min-hit-len 30] \
    [--merge-dist 100000] \
    [--min-array-len 10000] \
    [-h|--help]

Notes:
  - The actual output directory will be PREFIX_output_dir_base (e.g., AthCEN178_arrayout).
  - If --chr/--chr-file are omitted, the script auto-detects from split FASTAs
    (prefers chr1..chr5 if present, otherwise all sequences).

Examples:
  ./arrayfind.sh -g genome.fasta -c AthCEN178_consensus.fasta -o arrayout -p AthCEN178
  ./arrayfind.sh -g genome.fasta -c motif.fa -o results -p MyMotif --chr "chr1,chr3"
  ./arrayfind.sh -g genome.fasta -c motif.fa -o results -p MyMotif --chr-file chroms.txt
EOF
}

# ---- parse args ----
if [[ $# -eq 0 ]]; then print_help; exit 0; fi

GENOME=""
CONSENSUS=""
OUTBASE=""
PREFIX=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--genome) GENOME="${2:-}"; shift 2;;
    -c|--consensus) CONSENSUS="${2:-}"; shift 2;;
    -o|--outdir) OUTBASE="${2:-}"; shift 2;;
    -p|--prefix) PREFIX="${2:-}"; shift 2;;
    --chr) CHR_CSV="${2:-}"; shift 2;;
    --chr-file) CHR_FILE="${2:-}"; shift 2;;
    --min-hit-len) MIN_HIT_LEN="${2:-}"; shift 2;;
    --merge-dist) MERGE_DIST="${2:-}"; shift 2;;
    --min-array-len) MIN_ARRAY_LEN="${2:-}"; shift 2;;
    -h|--help) print_help; exit 0;;
    *) die "Unknown argument: $1 (use -h)";;
  esac
done

[[ -z "$GENOME"    ]] && { print_help; die "Provide -g/--genome"; }
[[ -z "$CONSENSUS" ]] && { print_help; die "Provide -c/--consensus"; }
[[ -z "$OUTBASE"   ]] && { print_help; die "Provide -o/--outdir"; }
[[ -z "$PREFIX"    ]] && { print_help; die "Provide -p/--prefix"; }

GENOME="$(readlink -f "$GENOME")"
CONSENSUS="$(readlink -f "$CONSENSUS")"

# ---- derive prefixed OUTDIR ----
OUTDIR="$(readlink -f "${OUTBASE}")"
mkdir -p "$OUTDIR"

# ---- load modules (optional) ----
ml SeqKit/2.9.0 || true
ml BLAST+/2.16.0-gompi-2024a || true
ml BEDTools/2.31.1-GCC-13.3.0 || true

# ---- sanity checks ----
need seqkit
need blastn
need bedtools
need awk
need cut
need basename
need sort
need xargs

# ---- output structure ----
SPLIT_DIR="$OUTDIR/split_seq"
BLAST_DIR="$OUTDIR/blast"
ARRAY_DIR="$OUTDIR"
mkdir -p "$SPLIT_DIR" "$BLAST_DIR" "$ARRAY_DIR"

# ---- split genome ----
if [[ ! -d "$SPLIT_DIR" || -z "$(ls -A "$SPLIT_DIR" 2>/dev/null || true)" ]]; then
  echo "[*] Splitting genome with seqkit..."
  seqkit split -i --by-id-prefix "" -O "$SPLIT_DIR" "$GENOME"
else
  echo "[*] Reusing existing split dir: $SPLIT_DIR"
fi

# seqkit sometimes creates a subfolder named split/
if [[ -d "$SPLIT_DIR/split" ]]; then
  SPLIT_DIR="$SPLIT_DIR/split"
fi

# ---- determine chromosome list ----
declare -a CHRS
if [[ -n "$CHR_FILE" ]]; then
  mapfile -t CHRS < <(grep -v '^\s*$' "$CHR_FILE")
elif [[ -n "$CHR_CSV" ]]; then
  IFS=',' read -r -a CHRS <<< "$CHR_CSV"
else
  if ls "$SPLIT_DIR"/chr{1,2,3,4,5}.fasta >/dev/null 2>&1; then
    CHRS=(chr1.fasta chr2.fasta chr3.fasta chr4.fasta chr5.fasta)
  else
    mapfile -t CHRS < <(ls "$SPLIT_DIR"/*.fasta | xargs -n1 basename)
  fi
fi

[[ ${#CHRS[@]} -eq 0 ]] && die "No chromosome FASTAs found/chosen."
echo "[*] Chromosomes to process: ${CHRS[*]}"

# ---- run BLAST and generate arrays ----
for chrfa in "${CHRS[@]}"; do
  chr_base="$(basename "$chrfa" .fasta)"
  chr_path="$SPLIT_DIR/$chrfa"
  [[ -f "$chr_path" ]] || die "Missing $chr_path"

  blast_txt="$BLAST_DIR/${chr_base}_${PREFIX}_blast.txt"
  merged_bed="$BLAST_DIR/${chr_base}_${PREFIX}_blast_filter_merge.bed"
  array_bed="$ARRAY_DIR/${chr_base}_${PREFIX}_array.bed"
  array_fa="$ARRAY_DIR/${chr_base}_${PREFIX}_array.fa"

  echo "[*] BLAST: $chr_base vs $CONSENSUS"
  blastn -subject "$chr_path" -query "$CONSENSUS" -outfmt 6 -max_target_seqs 5000000 > "$blast_txt"

  echo "[*] Filter+merge BLAST hits → BED ($chr_base)"
  # outfmt 6 columns: qseqid sseqid pident length qstart qend sstart send evalue bitscore
  awk -v MINLEN="$MIN_HIT_LEN" '$4>=MINLEN' "$blast_txt" \
    | cut -f2,9,10 \
    | awk '{if($3<$2){start=$3-1; end=$2}else{start=$2-1; end=$3}; print $1"\t"start"\t"end}' \
    | bedtools sort -i - \
    | bedtools merge -d "$MERGE_DIST" -i - \
    | awk -v MINARR="$MIN_ARRAY_LEN" '$3-$2>MINARR' \
    > "$merged_bed"

  cp "$merged_bed" "$array_bed"

  echo "[*] Extract array FASTA ($chr_base)"
  bedtools getfasta -fi "$chr_path" -bed "$merged_bed" -fo "$array_fa"
done

echo "[✓] Done."
echo "Output base (prefixed): $OUTDIR"
echo "  Split FASTAs:   $SPLIT_DIR"
echo "  BLAST results:  $BLAST_DIR/*_${PREFIX}_blast.txt"
echo "  Array BEDs:     $ARRAY_DIR/*_${PREFIX}_array.bed"
echo "  Array FASTAs:   $ARRAY_DIR/*_${PREFIX}_array.fa"