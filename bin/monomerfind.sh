#!/usr/bin/env bash
set -euo pipefail

# ======================================
# monomerfind.sh  (arrays-only)
# Run nhmmer on precomputed array FASTAs to extract monomers,
# filter by length, clean headers, and split into single-monomer FASTAs.
#
# Inputs expected in --arrays-dir:
#   <chr>_<PREFIX>_array.fa  (from arrayfind.sh)
#
# Outputs (under --outdir):
#   <chr>_<PREFIX>_array_hmmoutF.*                  (nhmmer outputs)
#   <chr>_<PREFIX>_array_monomer*.fa                (raw/filt/renamed)
#   <chr>_<PREFIX>_array_monomer_filt.format.fa.split/*.fa   (one file per monomer)
#   <PREFIX>_monomers/*.fa                          (all monomers, flat directory)
#
# Requires: nhmmer (HMMER), makehmmerdb, bedtools, bioawk, awk, sed, xargs
# ======================================

# ---------- defaults ----------
ARRAYS_DIR=""          # directory that contains *<chr>_<PREFIX>_array.fa
OUTDIR=""              # where to put monomer outputs (can be same as arrays dir)
PREFIX=""              # e.g., CEN178
HMMFILE=""             # profile HMM for nhmmer
MIN_MONOMER_LEN=120    # min monomer length to keep (bp)
CHR_CSV=""             # optional: comma-separated subset of chromosomes to process
CHR_FILE=""            # optional: file with one chromosome basename per line (e.g., chr1)

# ---------- helpers ----------
die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

print_help() {
  cat <<'EOF'
Usage:
  ./monomerfind.sh \
    --arrays-dir /path/to/arrays \
    --outdir /path/to/out \
    --prefix CEN178 \
    --hmm model.hmm \
    [--min-monomer-len 120] \
    [--chr "chr1,chr2,..."] \
    [--chr-file chr_list.txt]

Notes:
  - Expects array FASTAs named <chr>_<PREFIX>_array.fa inside --arrays-dir.
  - Outputs go under --outdir.
  - Adds two conveniences:
      1) <chr>_<PREFIX>_array_monomer_filt.format.fa.split/  (per-chr, one monomer per file)
      2) <PREFIX>_monomers/                                  (flat dir, one monomer per file)
EOF
}

# ---------- parse args ----------
[[ $# -eq 0 ]] && { print_help; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arrays-dir) ARRAYS_DIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --prefix) PREFIX="$2"; shift 2;;
    --hmm) HMMFILE="$2"; shift 2;;
    --min-monomer-len) MIN_MONOMER_LEN="$2"; shift 2;;
    --chr) CHR_CSV="$2"; shift 2;;
    --chr-file) CHR_FILE="$2"; shift 2;;
    -h|--help) print_help; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "$ARRAYS_DIR" ]] && die "Provide --arrays-dir"
[[ -z "$OUTDIR" ]]     && die "Provide --outdir"
[[ -z "$PREFIX" ]]     && die "Provide --prefix"
[[ -z "$HMMFILE" ]]    && die "Provide --hmm"

ARRAYS_DIR="$(readlink -f "$ARRAYS_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
HMMFILE="$(readlink -f "$HMMFILE")"
mkdir -p "$OUTDIR"

# ---- load modules (optional; ignore if not available) ----
if command -v ml >/dev/null 2>&1; then
  ml HMMER/3.4-gompi-2023a || true
  ml bioawk/1.0-GCC-12.3.0 || true
  ml BEDTools/2.31.1-GCC-13.3.0 || true
fi

# ---------- deps ----------
need nhmmer
need bedtools
need awk
need sed
need bioawk
need xargs
need makehmmerdb

# ---------- dirs ----------
MONO_DIR="$OUTDIR"
FLAT_SPLIT_DIR="$OUTDIR/${PREFIX}_monomers"   # <- flat directory with all single monomer FASTAs
mkdir -p "$MONO_DIR" "$FLAT_SPLIT_DIR"

# ---------- determine which array FASTAs to use ----------
declare -a CHRS
if [[ -n "$CHR_FILE" ]]; then
  mapfile -t CHRS < <(grep -v '^\s*$' "$CHR_FILE")
elif [[ -n "$CHR_CSV" ]]; then
  IFS=',' read -r -a CHRS <<< "$CHR_CSV"
else
  # infer from files present
  mapfile -t CHRS < <(ls "$ARRAYS_DIR"/*_"$PREFIX"_array.fa 2>/dev/null \
                       | xargs -n1 basename \
                       | sed -E "s/_${PREFIX}_array\.fa$//" )
fi
[[ ${#CHRS[@]} -eq 0 ]] && die "No array FASTAs found for prefix '$PREFIX' in $ARRAYS_DIR"

echo "[*] Will process arrays for: ${CHRS[*]}"

# ---------- nhmmer over arrays ----------
for chr in "${CHRS[@]}"; do
  array_fa="$ARRAYS_DIR/${chr}_${PREFIX}_array.fa"
  [[ -f "$array_fa" ]] || { echo "[!] Skip $chr (missing $array_fa)"; continue; }
  if [[ ! -s "$array_fa" ]]; then
    echo "[!] Skip $chr (empty array fasta)."
    continue
  fi

  echo "[*] nhmmer: $chr"
  db="$ARRAYS_DIR/${chr}_${PREFIX}_array.fa.db"
  makehmmerdb "$array_fa" "$db"

  tbl="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.tbl"
  out="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.out"
  nhmmer --tblout "$tbl" -o "$out" "$HMMFILE" "$db"

  # Parse coordinates from the human-readable output
  bed="$MONO_DIR/${chr}_${PREFIX}_array_hmmoutF.out.bed"
  awk '/Scores for complete hits:/, /^$/' "$out" \
    | awk '{ gsub(/^    /,"",$0); print }' \
    | awk '$1 ~ /^[0-9]/ {
        chr=$4; start=$5-1; end=$6;
        if (end < start) { t=start; start=end; end=t }
        print chr"\t"start"\t"end
      }' \
    | bedtools sort -i - \
    > "$bed"

  # Extract monomers, filter by length, clean headers
  mono_raw="$MONO_DIR/${chr}_${PREFIX}_array_monomer.fa"
  bedtools getfasta -fi "$array_fa" -bed "$bed" -fo "$mono_raw"

  mono_filt="$MONO_DIR/${chr}_${PREFIX}_array_monomer_filt.fa"
  bioawk -c fastx -v L="$MIN_MONOMER_LEN" 'length($seq) >= L { print ">"$name; print $seq }' \
    "$mono_raw" > "$mono_filt"

  # Normalize headers: replace ':' and '-' with '_' so they’re filename-safe
  mono_fmt="$MONO_DIR/${chr}_${PREFIX}_array_monomer_filt.format.fa"
  sed 's/:/_/g; s/-/_/g' "$mono_filt" > "$mono_fmt"

  # ---- NEW: split into one FASTA per monomer ----
  #   A) per-chromosome split folder (compat with other tools)
  #   B) flat folder OUTDIR/${PREFIX}_monomers with all monomers
  split_dir="$MONO_DIR/${chr}_${PREFIX}_array_monomer_filt.format.fa.split"
  mkdir -p "$split_dir"

  bioawk -c fastx -v D1="$split_dir" -v D2="$FLAT_SPLIT_DIR" '
    {
      # sanitize header further for filenames (just in case)
      nm=$name
      gsub(/[^A-Za-z0-9_.-]/,"_",nm)
      fn = nm ".fa"
      print ">"$name"\n"$seq > (D1 "/" fn)
      print ">"$name"\n"$seq > (D2 "/" fn)
    }
  ' "$mono_fmt"
done

echo "[✓] Done."
echo "Monomer outputs:           $MONO_DIR"
echo "Per-monomer (per-chr):     <chr>_${PREFIX}_array_monomer_filt.format.fa.split/*.fa"
echo "Per-monomer (flat):        $FLAT_SPLIT_DIR/*.fa" 