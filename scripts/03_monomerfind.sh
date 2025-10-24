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
# HiReNET: monomerfind.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Detect and extract individual monomers within tandem repeat
#   arrays identified by arrayfind. Runs nhmmer using the provided
#   profile HMM, converts array-relative hits to chromosome 
#   coordinates, and outputs per-chromosome monomer FASTAs.
#
# Dependencies:
#   - nhmmer (HMMER)
#   - seqkit
#   - bedtools
#
# Input:
#   - Array FASTA files from HiReNET arrayfind output
#   - Chromosome FASTAs for coordinate mapping
#   - Profile HMM model (from HiReNET getphmm)
#
# Output:
#   - Monomer hit tables (per chromosome)
#   - Extracted monomer FASTAs (per chromosome)
#
# Example:
#   HiReNET monomerfind \
#       --arrays-dir AthCEN178_arrayout \
#       --chrom-dir  AthCEN178_arrayout/split_seq \
#       --outdir     AthCEN178_monomerout \
#       --prefix     AthCEN178 \
#       --hmm        AthCEN178_phmm/AthCEN178.hmm \
#       --chr "chr1,chr2,chr3,chr4,chr5"
#
# Notes:
#   - Automatically converts relative monomer coordinates to genomic
#     positions using array headers (chr:start-end format).
#   - Works with single or multiple chromosomes.
# ==============================================================

# ---------- defaults ----------
ARRAYS_DIR=""
CHROM_DIR=""
OUTDIR=""
PREFIX=""
HMMFILE=""
MIN_MONOMER_LEN=120
CHR_CSV=""
CHR_FILE=""
# ---------- helpers ----------
die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

# ---- usage / help ----
usage() {
  cat <<'EOF'
============================================================
                     HiReNET: monomerfind
============================================================

Description:
  Detect individual monomers within tandem repeat arrays by 
  scanning with a profile HMM. Converts array-relative hits 
  into chromosome coordinates and extracts corresponding sequences.

------------------------------------------------------------
Usage:
  HiReNET monomerfind \
      --arrays-dir <path/to/arrays> \
      --chrom-dir  <path/to/chrom_fastas> \
      --outdir     <path/to/output_dir> \
      --prefix     <prefix> \
      --hmm        <model.hmm> [options]

Required arguments:
  --arrays-dir <dir>       Directory of array FASTAs (from arrayfind)
  --chrom-dir  <dir>       Directory of chromosome FASTAs
  --outdir     <dir>       Output directory (created if missing)
  --prefix     <name>      Prefix for output files
  --hmm        <file>      Profile HMM (from getphmm output)

Optional arguments:
  --min-monomer-len <int>  Minimum monomer length to keep (default: 120)
  --chr "chr1,chr2,..."    Comma-separated chromosome names
  --chr-file <file>        File containing chromosome list (one per line)
  -h, --help               Show this help message and exit

------------------------------------------------------------
Notes:
  • nhmmer runs against array FASTAs to detect repeat monomers.
  • Coordinates are lifted to genome scale using array headers.
  • Monomer FASTAs are extracted per chromosome using seqkit.

Example:
  HiReNET monomerfind \
      --arrays-dir AthCEN178_arrayout \
      --chrom-dir  AthCEN178_arrayout/split_seq \
      --outdir     AthCEN178_monomerout \
      --prefix     AthCEN178 \
      --hmm        AthCEN178_phmm/AthCEN178.hmm \
      --chr "chr1,chr2,chr3,chr4,chr5"

============================================================
EOF
}

# Backward compatibility (in case older code still calls print_help)
print_help() { usage; }

# ---------- parse args ----------
[[ $# -eq 0 ]] && { print_help; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arrays-dir) ARRAYS_DIR="$2"; shift 2;;
    --chrom-dir)  CHROM_DIR="$2"; shift 2;;
    --outdir)     OUTDIR="$2"; shift 2;;
    --prefix)     PREFIX="$2"; shift 2;;
    --hmm)        HMMFILE="$2"; shift 2;;
    --min-monomer-len) MIN_MONOMER_LEN="$2"; shift 2;;
    --chr)        CHR_CSV="$2"; shift 2;;
    --chr-file)   CHR_FILE="$2"; shift 2;;
    -h|--help)    print_help; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "$ARRAYS_DIR" ]] && die "Provide --arrays-dir"
[[ -z "$CHROM_DIR"  ]] && die "Provide --chrom-dir (folder with chr FASTAs)"
[[ -z "$OUTDIR"     ]] && die "Provide --outdir"
[[ -z "$PREFIX"     ]] && die "Provide --prefix"
[[ -z "$HMMFILE"    ]] && die "Provide --hmm"

ARRAYS_DIR="$(readlink -f "$ARRAYS_DIR")"
CHROM_DIR="$(readlink -f "$CHROM_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
HMMFILE="$(readlink -f "$HMMFILE")"
mkdir -p "$OUTDIR"

# ---------- deps (Linux PATH-based) ----------
need nhmmer
need makehmmerdb
need bedtools
need awk
need sed
need bioawk
need xargs

# ---------- dirs ----------
MONO_DIR="$OUTDIR"
FLAT_SPLIT_DIR="$OUTDIR/${PREFIX}_monomers"   # flat directory with all single monomer FASTAs (genomic)
mkdir -p "$MONO_DIR" "$FLAT_SPLIT_DIR"

# ---------- determine chromosomes from arrays ----------
declare -a CHRS
if [[ -n "$CHR_FILE" ]]; then
  mapfile -t CHRS < <(grep -v '^[[:space:]]*$' "$CHR_FILE")
elif [[ -n "$CHR_CSV" ]]; then
  IFS=',' read -r -a CHRS <<< "$CHR_CSV"
else
  mapfile -t CHRS < <(ls "$ARRAYS_DIR"/*_"$PREFIX"_array.fa 2>/dev/null \
                      | xargs -n1 basename \
                      | sed -E "s/_${PREFIX}_array\.fa$//" )
fi
[[ ${#CHRS[@]} -eq 0 ]] && die "No array FASTAs found for prefix '$PREFIX' in $ARRAYS_DIR"

echo "[*] Arrays for: ${CHRS[*]}"

# ---------- helper: find chr fasta path ----------
chr_fa_path() {
  local c="$1"
  if   [[ -f "$CHROM_DIR/${c}.fasta" ]]; then echo "$CHROM_DIR/${c}.fasta"
  elif [[ -f "$CHROM_DIR/${c}.fa"    ]]; then echo "$CHROM_DIR/${c}.fa"
  else echo ""; fi
}

# ---------- nhmmer over arrays, then lift to chr coords, then getfasta from chr ----------
for chr in "${CHRS[@]}"; do
  array_fa="$ARRAYS_DIR/${chr}_${PREFIX}_array.fa"
  [[ -f "$array_fa" ]] || { echo "[!] Skip $chr (missing $array_fa)"; continue; }
  [[ -s "$array_fa" ]] || { echo "[!] Skip $chr (empty array fasta)."; continue; }

  chr_fa="$(chr_fa_path "$chr")"
  [[ -n "$chr_fa" ]] || { echo "[!] Skip $chr (missing ${chr}.fa|fasta in $CHROM_DIR)"; continue; }

  echo "[*] nhmmer (array): $chr"
  # run nhmmer directly on the array fasta (hmmer DB improves speed on large inputs)
  tbl="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.tbl"
  out="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.out"
  db="$ARRAYS_DIR/${chr}_${PREFIX}_array.fa.db"
  makehmmerdb "$array_fa" "$db"
  nhmmer --tblout "$tbl" -o "$out" "$HMMFILE" "$db"

  # Parse array-relative coords from the human-readable output
  arr_bed="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.out.bed"
  awk '/Scores for complete hits:/, /^$/' "$out" \
    | awk '{ gsub(/^    /,"",$0); print }' \
    | awk '$1 ~ /^[0-9]/ {
        arr=$4;
        start=$5-1; end=$6;   # array-relative, convert to 0-based half-open
        if (end < start) { t=start; start=end; end=t }
        print arr"\t"start"\t"end
      }' > "$arr_bed"

  # Convert array-relative to chromosome coords:
  gen_bed="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.out.genomic.bed"
  awk 'BEGIN{FS=OFS="\t"}
       {
         split($1, A, ":"); chr=A[1];
         split(A[2], R, "-"); AS=R[1]+0; AE=R[2]+0;
         relS=$2+0; relE=$3+0;
         gS = AS + relS; gE = AS + relE;
         print chr, gS, gE
       }' "$arr_bed" \
    | bedtools sort -i - > "$gen_bed"

  gen_bed_filt="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.out.genomic.filt${MIN_MONOMER_LEN}.bed"
  awk -v L="$MIN_MONOMER_LEN" 'BEGIN{FS=OFS="\t"} ($3-$2) >= L' "$gen_bed" \
    | bedtools sort -i - > "$gen_bed_filt"

  # ---- Keep array-based monomers (optional, unchanged from original) ----
  mono_raw_array="$MONO_DIR/${chr}_${PREFIX}_array_monomer.fa"
  bedtools getfasta -fi "$array_fa" -bed "$arr_bed" -fo "$mono_raw_array" 2>/dev/null || true

  mono_filt_array="$MONO_DIR/${chr}_${PREFIX}_array_monomer_filt.fa"
  if [[ -s "$mono_raw_array" ]]; then
    bioawk -c fastx -v L="$MIN_MONOMER_LEN" 'length($seq) >= L { print ">"$name; print $seq }' \
      "$mono_raw_array" > "$mono_filt_array" || true

    mono_fmt_array="$MONO_DIR/${chr}_${PREFIX}_array_monomer_filt.format.fa"
    sed 's/:/_/g; s/-/_/g' "$mono_filt_array" > "$mono_fmt_array" || true

    split_dir_array="$MONO_DIR/${chr}_${PREFIX}_array_monomer_filt.format.fa.split"
    mkdir -p "$split_dir_array"

    [[ -s "$mono_fmt_array" ]] && bioawk -c fastx -v D="$split_dir_array" '
      { nm=$name; gsub(/[^A-Za-z0-9_.-]/,"_",nm);
        print ">"$name"\n"$seq > (D "/" nm ".fa")
      }' "$mono_fmt_array" || true
  fi

  # NOW: extract monomers from chromosome FASTA using genomic BED
  mono_gen="$MONO_DIR/${chr}_${PREFIX}_monomer_genomic.fa"
  bedtools getfasta -fi "$chr_fa" -bed "$gen_bed_filt" -fo "$mono_gen"

  # Keep short header style; just normalize ':' and '-' to '_'
  mono_fmt_gen="$MONO_DIR/${chr}_${PREFIX}_monomer_genomic_filt.format.fa"
  sed 's/:/_/g; s/-/_/g' "$mono_gen" > "$mono_fmt_gen"

  # Split per monomer (short names) into per-chr folder and flat folder
  split_dir_gen="$MONO_DIR/${chr}_${PREFIX}_monomer_genomic_filt.format.fa.split"
  mkdir -p "$split_dir_gen" "$FLAT_SPLIT_DIR"

  # Stream-safe splitter (closes files) – keeps short names
  awk -v D1="$split_dir_gen" -v D2="$FLAT_SPLIT_DIR" '
    BEGIN{ RS=">"; ORS=""; }
    NR>1{
      n = split($0, L, /\r?\n/);
      hdr = L[1];
      gsub(/[^A-Za-z0-9_.:-]/,"_",hdr);   # keep short form like chr1:123-456 → chr1_123_456 after sed above
      seq = "";
      for (i=2; i<=n; i++) seq = seq L[i];
      if (length(seq)==0) next;

      # per-chr split
      fn1 = D1 "/" hdr ".fa";
      print ">" hdr "\n" seq "\n" > fn1; close(fn1);

      # flat dir (same short filename)
      fn2 = D2 "/" hdr ".fa";
      print ">" hdr "\n" seq "\n" > fn2; close(fn2);
    }
  ' "$mono_fmt_gen"
done

echo "[✓] Done."
echo "Arrays dir:                $ARRAYS_DIR"
echo "Chrom FASTAs dir:          $CHROM_DIR"
echo "Array-level nhmmer:        <chr>_${PREFIX}_array_hmmoutF.tbl/.out"
echo "Genomic BED (lifted):      <chr>_${PREFIX}_array_hmmoutF.out.genomic.bed"
echo "Genomic monomers (FA):     <chr>_${PREFIX}_monomer_genomic*.fa"
echo "Per-monomer (flat):        $FLAT_SPLIT_DIR/*.fa"
