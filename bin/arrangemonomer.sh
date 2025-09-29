#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# arrangemonomer.sh (pre-split monomers ONLY; flat folder)
#
# Input:
#   --monomer-dir  /path/<PREFIX>_monomers/*.fa
#     Each file is a single monomer, named like:
#       chr1_14389697_18147662_5086_5263.fa
#       ^chr ^array_start ^array_end  ^mon_start ^mon_end
#
# Outputs (under --outdir):
#   <chr>_bins_<BIN>bp_bed/*.bed
#   <chr>_<PREFIX>_array_monomer/*.fa                 (per-bin FASTAs)
#   <chr>_<PREFIX>_array_monomer_filt.<BIN>bp_with_bins.bed
#   ${PREFIX}_monomer_bed_inbin.txt                   (combined from all chr)
#   ${PREFIX}_bin_monomers/                           (all per-bin FASTAs)
#
# Usage:
#   ./arrangemonomer.sh \
#     --monomer-dir /scratch/.../AthCEN178_monomers \
#     --outdir      /scratch/.../arrange_out \
#     --prefix      AthCEN178 \
#     [--bin 10000] \
#     [--chr "chr1,chr2,..."] | [--chr-file list.txt]
#
# Requires: awk, sort, sed, find, xargs, mkdir, basename
# =========================================================

# ---------- args ----------
MONOMER_DIR=""
OUTDIR=""
PREFIX=""
BIN=10000
CHR_CSV=""
CHR_FILE=""

usage() {
  cat <<USAGE
Usage:
  $0 --monomer-dir DIR --outdir DIR --prefix NAME [--bin N] [--chr "chr1,chr2"] [--chr-file file.txt]

Notes:
  - This script expects a flat directory of single-monomer FASTAs:
      <MONOMER_DIR>/*.fa
    with filenames like 'chr1_14389697_18147662_5086_5263.fa'.
USAGE
  exit 1
}

[[ $# -eq 0 ]] && usage
while [[ $# -gt 0 ]]; do
  case "$1" in
    --monomer-dir) MONOMER_DIR="${2:-}"; shift 2;;
    --outdir)      OUTDIR="${2:-}"; shift 2;;
    --prefix)      PREFIX="${2:-}"; shift 2;;
    --bin)         BIN="${2:-}"; shift 2;;
    --chr)         CHR_CSV="${2:-}"; shift 2;;
    --chr-file)    CHR_FILE="${2:-}"; shift 2;;
    -h|--help)     usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

# ---------- validate ----------
die() { echo "ERROR: $*" >&2; exit 2; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing <$1>"; }
[[ -z "$MONOMER_DIR" ]] && die "--monomer-dir required"
[[ -z "$OUTDIR"      ]] && die "--outdir required"
[[ -z "$PREFIX"      ]] && die "--prefix required"

need awk; need sort; need sed; need find; need xargs; need basename; need mkdir

MONOMER_DIR="$(readlink -f "$MONOMER_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"

# ---------- choose chromosomes ----------
declare -a CHRS
if [[ -n "$CHR_FILE" ]]; then
  mapfile -t CHRS < <(grep -v '^\s*$' "$CHR_FILE")
elif [[ -n "$CHR_CSV" ]]; then
  IFS=',' read -r -a CHRS <<< "$CHR_CSV"
else
  # infer from filenames
  mapfile -t CHRS < <(find "$MONOMER_DIR" -maxdepth 1 -type f -name '*.fa' -printf '%f\n' \
                      | awk -F'_' '{print $1}' | sort -u)
fi
[[ ${#CHRS[@]} -eq 0 ]] && die "No chromosomes inferred/found."

echo "[*] Monomer dir : $MONOMER_DIR"
echo "[*] Outdir      : $OUTDIR"
echo "[*] Prefix      : $PREFIX"
echo "[*] BIN         : $BIN"
echo "[*] Chr list    : ${CHRS[*]}"

# ---------- per-chromosome processing ----------
COMBINED_ALL="$OUTDIR/${PREFIX}_monomer_bed_inbin.txt"
: > "$COMBINED_ALL"

for chr in "${CHRS[@]}"; do
  echo "[*] Processing $chr"

  # 1) Build a 5-col "bed-like" table by parsing filenames in the monomer dir.
  #    Expected filename: chr_arrayStart_arrayEnd_monoStart_monoEnd.fa
  bed5="$OUTDIR/${chr}_${PREFIX}_array_monomer_filt.format.bed"
  : > "$bed5"

  # collect files for this chr
  mapfile -t chr_files < <(find "$MONOMER_DIR" -maxdepth 1 -type f -name "${chr}_*.fa" -printf '%f\n' | sort)
  if [[ ${#chr_files[@]} -eq 0 ]]; then
    echo "  [!] No monomer files for $chr under $MONOMER_DIR — skipping."
    continue
  fi

  # Parse into 5 cols
  #   chr array_start array_end monomer_start monomer_end
  for f in "${chr_files[@]}"; do
    base="${f%.fa}"
    # split on "_"
    IFS='_' read -r c astart aend mstart mend <<< "$base" || true
    if [[ -z "${c:-}" || -z "${astart:-}" || -z "${aend:-}" || -z "${mstart:-}" || -z "${mend:-}" ]]; then
      echo "  [!] Skip malformed filename: $f"
      continue
    fi
    echo -e "${c}\t${astart}\t${aend}\t${mstart}\t${mend}" >> "$bed5"
  done

  # If nothing parsed, skip chr
  if [[ ! -s "$bed5" ]]; then
    echo "  [!] No usable monomer entries for $chr — skipping."
    continue
  fi

  # 2) Create per-bin BEDs (array-relative binning)
  OUTDIR_BED="$OUTDIR/${chr}_bins_${BIN}bp_bed"
  mkdir -p "$OUTDIR_BED"

  COMBINED_CHR="$OUTDIR/${chr}_${PREFIX}_array_monomer_filt.${BIN}bp_with_bins.bed"
  : > "$COMBINED_CHR"

  awk -v BIN="$BIN" -v OUTDIR="$OUTDIR_BED" -v COMB="$COMBINED_CHR" 'BEGIN{OFS="\t"}
    {
      chr=$1; astart=$2+0; aend=$3+0; mstart=$4+0; mend=$5+0
      arraysz = aend - astart; if (arraysz <= 0) next

      # array-relative binning (monomer coords are already array coords if produced by your pipeline)
      bstart = int(mstart / BIN) * BIN
      bend   = bstart + BIN
      if (bend > arraysz) bend = arraysz

      fname = OUTDIR "/" chr "_" astart "_" aend "_" bstart "_" bend ".bed"
      print chr, astart, aend, mstart, mend >> fname
      print chr, astart, aend, mstart, mend, bstart, bend >> COMB
    }
  ' "$bed5"

  # Sort per-bin BEDs and combined file
  for f in "$OUTDIR_BED"/*.bed; do
    [[ -e "$f" ]] || continue
    sort -k1,1 -k2,2n -k4,4n "$f" | uniq > "${f}.sorted" && mv "${f}.sorted" "$f"
  done
  sort -k1,1 -k2,2n -k6,6n "$COMBINED_CHR" -o "$COMBINED_CHR"

  # 3) Build per-bin FASTAs by concatenating the single-monomer files
  BIN_FASTA_DIR="$OUTDIR/${chr}_${PREFIX}_array_monomer"
  mkdir -p "$BIN_FASTA_DIR"

  while IFS= read -r -d '' bed; do
    base="$(basename "$bed" .bed)"
    outfa="$BIN_FASTA_DIR/${base}.fa"

    # monomer IDs (filenames without .fa) are exactly the 5 fields joined by "_"
    # Build paths and concat
    awk 'BEGIN{OFS="_"} {print $1,$2,$3,$4,$5}' "$bed" \
      | while read -r id; do
          mf="$MONOMER_DIR/${id}.fa"
          [[ -f "$mf" ]] && cat "$mf"
        done > "$outfa"
  done < <(find "$OUTDIR_BED" -type f -name '*.bed' -print0)

  # Append to global combined
  cat "$COMBINED_CHR" >> "$COMBINED_ALL"
done


# 4) Collect all per-bin FASTAs into ${PREFIX}_bin_monomers/
BIN_FLAT_DIR="$OUTDIR/${PREFIX}_bin_monomers"
mkdir -p "$BIN_FLAT_DIR"

# Copy every per-bin FASTA produced under .../<chr>_${PREFIX}_array_monomer/*.fa
find "$OUTDIR" -type f -path "*/*_${PREFIX}_array_monomer/*.fa" -print0 \
  | xargs -0 -I{} cp -f "{}" "$BIN_FLAT_DIR/"

  
echo "[✓] Done."
echo "  Combined bin table : $COMBINED_ALL"
echo "  Per-bin FASTAs     : <chr>_${PREFIX}_array_monomer/*.fa"
echo "  All bin FASTAs     : $BIN_FLAT_DIR/"