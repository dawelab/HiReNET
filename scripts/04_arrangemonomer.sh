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
# HiReNET: arrangemonomer.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Organize monomer coordinates into array-anchored genomic bins
#   of a specified size (e.g., 10 kb). Combines array, monomer, and
#   BED-based coordinate data to produce bin-level monomer summaries
#   and per-bin FASTAs for downstream comparison and classification.
#
# Dependencies:
#   - bedtools
#   - awk
#   - seqkit
#
# Input:
#   --arrays-dir        Directory with per-chromosome array BEDs  
#                       (e.g., chrX_<PREFIX>_array.bed)
#   --genomic-bed-dir   Directory with genomic monomer BEDs  
#                       (e.g., chrX_<PREFIX>_array_hmmoutF.out.genomic.bed)
#   --monomer-dir       Directory with single-monomer FASTAs  
#                       (e.g., chr_arrayStart_arrayEnd_monoStart_monoEnd.fa)
#
# Output (under --outdir):
#   - <chr>_bins_<BIN>bp_bed/*.bed
#   - <chr>_<PREFIX>_array_monomer/*.fa
#   - <chr>_<PREFIX>_array_monomer_filt.<BIN>bp_with_bins.bed
#   - ${PREFIX}_monomer_bed_inbin.txt
#   - ${PREFIX}_bin_monomers/
#
# Example:
#   HiReNET arrangemonomer \
#       --arrays-dir AthCEN178_arrayout \
#       --genomic-bed-dir AthCEN178_monomerout \
#       --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
#       --outdir AthCEN178_arrangemonomer_10kb \
#       --prefix AthCEN178 \
#       --bin 10000 \
#       --chr "chr1,chr2,chr3,chr4,chr5"
#
# Notes:
#   - Automatically builds per-bin BEDs and FASTAs for each chromosome.
#   - The --bin argument controls bin length (default: 10000 bp).
#   - Requires consistent file naming from previous HiReNET stages.
# ==============================================================

# ---------- args ----------
ARRAYS_DIR=""
GENBED_DIR=""
MONOMER_DIR=""
OUTDIR=""
PREFIX=""
BIN=10000
MINCNT=5
CHR_CSV=""
CHR_FILE=""

usage() {
  cat <<'EOF'
============================================================
                     HiReNET: arrangemonomer
============================================================

Description:
  Build array-anchored, fixed-size genomic bins from monomer BEDs
  and FASTAs. Combines array, monomer, and BED coordinate data to
  produce per-bin monomer summaries and sequences for downstream
  comparison and classification.

------------------------------------------------------------
Usage:
  HiReNET arrangemonomer \
      --arrays-dir <path/to/arrays> \
      --genomic-bed-dir <path/to/genomic_beds> \
      --monomer-dir <path/to/monomers> \
      --outdir <path/to/output_dir> \
      --prefix <name> [options]

Required arguments:
  --arrays-dir <dir>        Directory with per-chromosome array BEDs
  --genomic-bed-dir <dir>   Directory with genomic monomer BEDs
  --monomer-dir <dir>       Directory with monomer FASTAs
  --outdir <dir>            Output directory (created if missing)
  --prefix <name>           Prefix for output files

Optional arguments:
  --bin <int>               Bin size in base pairs (default: 10000)
  --mincnt <int>            Minimum monomers per bin (default: 5)
  --chr "chr1,chr2"         Comma-separated list of chromosomes
  --chr-file <file>         File containing chromosome list (one per line)
  -h, --help                Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • <chr>_bins_<BIN>bp_bed/*.bed
  • <chr>_<PREFIX>_array_monomer/*.fa
  • <chr>_<PREFIX>_array_monomer_filt.<BIN>bp_with_bins.bed
  • ${PREFIX}_monomer_bed_inbin.txt
  • ${PREFIX}_bin_monomers/

Example:
  HiReNET arrangemonomer \
      --arrays-dir AthCEN178_arrayout \
      --genomic-bed-dir AthCEN178_monomerout \
      --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
      --outdir AthCEN178_arrangemonomer_10kb \
      --prefix AthCEN178 \
      --bin 10000 \
      --chr "chr1,chr2,chr3,chr4,chr5"

============================================================
EOF
}

# Backward compatibility (if older calls use print_help)
print_help() { usage; }

[[ $# -eq 0 ]] && { usage; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --arrays-dir)       ARRAYS_DIR="${2:-}"; shift 2;;
    --genomic-bed-dir)  GENBED_DIR="${2:-}"; shift 2;;
    --monomer-dir)      MONOMER_DIR="${2:-}"; shift 2;;
    --outdir)           OUTDIR="${2:-}"; shift 2;;
    --prefix)           PREFIX="${2:-}"; shift 2;;
    --bin)              BIN="${2:-}"; shift 2;;
    --mincnt)           MINCNT="${2:-}"; shift 2;;
    --chr)              CHR_CSV="${2:-}"; shift 2;;
    --chr-file)         CHR_FILE="${2:-}"; shift 2;;
    -h|--help)          usage; exit 0;;
    *) echo "Unknown argument: $1"; usage; exit 1;;
  esac
done

# ---------- validate ----------
die()  { echo "ERROR: $*" >&2; exit 2; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing <$1> (activate env?)"; }

[[ -z "$ARRAYS_DIR"  ]] && die "--arrays-dir required"
[[ -z "$GENBED_DIR"  ]] && die "--genomic-bed-dir required"
[[ -z "$MONOMER_DIR" ]] && die "--monomer-dir required"
[[ -z "$OUTDIR"      ]] && die "--outdir required"
[[ -z "$PREFIX"      ]] && die "--prefix required"

# Dependencies (Linux PATH-based)
need awk; need sort; need bedtools; need find; need xargs; need basename; need mkdir

# Absolute paths (Linux)
ARRAYS_DIR="$(readlink -f "$ARRAYS_DIR")"
GENBED_DIR="$(readlink -f "$GENBED_DIR")"
MONOMER_DIR="$(readlink -f "$MONOMER_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"

# ---------- chromosome list ----------
declare -a CHRS
if [[ -n "$CHR_FILE" ]]; then
  mapfile -t CHRS < <(grep -v '^[[:space:]]*$' "$CHR_FILE")
elif [[ -n "$CHR_CSV" ]]; then
  IFS=',' read -r -a CHRS <<< "$CHR_CSV"
else
  # infer chr names from array BEDs
  mapfile -t CHRS < <(find "$ARRAYS_DIR" -maxdepth 1 -type f -name "chr*_${PREFIX}_array.bed" -printf '%f\n' \
                      | sed -E "s/_${PREFIX}_array\.bed$//" | sort -u)
fi
[[ ${#CHRS[@]} -eq 0 ]] && die "No chromosomes inferred (check --arrays-dir and --prefix)."

echo "[*] Arrays dir     : $ARRAYS_DIR"
echo "[*] Genomic BED dir: $GENBED_DIR"
echo "[*] Monomer dir    : $MONOMER_DIR"
echo "[*] Outdir         : $OUTDIR"
echo "[*] Prefix         : $PREFIX"
echo "[*] BIN            : $BIN"
echo "[*] MINCNT         : $MINCNT"
echo "[*] Chr list       : ${CHRS[*]}"

# ---------- combined output ----------
COMBINED_ALL="$OUTDIR/${PREFIX}_monomer_bed_inbin.txt"
: > "$COMBINED_ALL"

# ---------- per-chromosome ----------
for chr in "${CHRS[@]}"; do
  echo "[*] Processing $chr"

  arrays_bed="$ARRAYS_DIR/${chr}_${PREFIX}_array.bed"
  genbed="$GENBED_DIR/${chr}_${PREFIX}_array_hmmoutF.out.genomic.filt120.bed"

  [[ -s "$arrays_bed" ]] || { echo "  [!] Missing/empty $arrays_bed — skip $chr"; continue; }
  [[ -s "$genbed"    ]] || { echo "  [!] Missing/empty $genbed — skip $chr"; continue; }

  # 1) Make bed5: chr  array_start  array_end  monomer_start(genomic)  monomer_end(genomic)
  bed5="$OUTDIR/${chr}_${PREFIX}_array_monomer_filt.format.bed"
  : > "$bed5"

  bedtools intersect -a "$genbed" -b "$arrays_bed" -wa -wb \
    | awk -v OFS='\t' '
        {
          chrA=$1; mS=$2+0; mE=$3+0;
          chrB=$4; aS=$5+0; aE=$6+0;
          if (chrA!=chrB) next;
          # keep monomers fully inside the array interval
          if (mS < aS || mE > aE) next;
          print chrA, aS, aE, mS, mE;
        }
      ' \
    | sort -k1,1 -k2,2n -k3,3n -k4,4n \
    > "$bed5"

  [[ -s "$bed5" ]] || { echo "  [!] No monomers within arrays for $chr"; continue; }

  # 2) Bin (genomic coords; anchored at first monomer’s genomic start per array)
  OUTDIR_BED="$OUTDIR/${chr}_bins_${BIN}bp_bed"
  mkdir -p "$OUTDIR_BED"

  COMBINED_CHR="$OUTDIR/${chr}_${PREFIX}_array_monomer_filt.${BIN}bp_with_bins.bed"
  : > "$COMBINED_CHR"

  awk -v BIN="$BIN" -v MINCNT="$MINCNT" -v OUTDIR="$OUTDIR_BED" -v COMB="$COMBINED_CHR" '
    BEGIN{FS=OFS="\t"}
    {
      chr=$1; aS=$2+0; aE=$3+0; gS=$4+0; gE=$5+0
      key = chr "|" aS "|" aE
      n[key]++
      GS[key, n[key]] = gS
      REC[key, n[key]] = $0
      if (!(key in lo) || gS < lo[key]) lo[key] = gS
      if (!(key in hi) || gS > hi[key]) hi[key] = gS
    }
    function sort_idx_by_GS(key, idx, tmp, i, N) {
      delete tmp; delete idx
      for (i=1; i<=n[key]; i++) tmp[i] = GS[key,i]
      N = asorti(tmp, idx, "@val_num_asc")
      return N
    }
    END{
      for (key in n){
        split(key, K, "|"); chr=K[1]; aS=K[2]; aE=K[3]
        N = sort_idx_by_GS(key, IDX); if (N==0) continue

        start = lo[key]; stop = hi[key]

        for (b = start; b <= stop + BIN; b += BIN) {
          bS = b; bE = b + BIN
          cnt = 0
          for (i=1; i<=N; i++) {
            gs = GS[key, IDX[i]]
            if (gs >= bS && gs < bE) cnt++
          }
          if (cnt > MINCNT) {
            fname = OUTDIR "/" chr "_" aS "_" aE "_" bS "_" bE ".bed"
            for (i=1; i<=N; i++) {
              rec = REC[key, IDX[i]]
              split(rec, F, FS)  # F: chr aS aE gS gE
              gs = F[4]+0
              if (gs >= bS && gs < bE) {
                print rec >> fname
                print rec, bS, bE >> COMB
              }
            }
          }
        }
      }
    }
  ' "$bed5"

  # tidy
  for f in "$OUTDIR_BED"/*.bed; do
    [[ -e "$f" ]] || continue
    sort -k1,1 -k2,2n -k4,4n "$f" | uniq > "${f}.sorted" && mv "${f}.sorted" "$f"
  done
  sort -k1,1 -k2,2n -k6,6n "$COMBINED_CHR" -o "$COMBINED_CHR"
  cat "$COMBINED_CHR" >> "$COMBINED_ALL"
done

# 3) Build per-bin FASTAs for every chromosome, then gather into one flat folder
for chr in "${CHRS[@]}"; do
  echo "[*] Processing $chr monomers"
  OUTDIR_BED="$OUTDIR/${chr}_bins_${BIN}bp_bed"
  [ -d "$OUTDIR_BED" ] || continue

  BIN_FASTA_DIR="$OUTDIR/${chr}_${PREFIX}_array_monomer"
  mkdir -p "$BIN_FASTA_DIR"

  for bed in "$OUTDIR_BED"/*.bed; do
    [ -e "$bed" ] || continue
    base="$(basename "$bed" .bed)"
    outfa="$BIN_FASTA_DIR/${base}.fa"

    # bed columns: chr aStart aEnd gStart gEnd
    awk -v MONO="$MONOMER_DIR" '
      {
        chr=$1; aS=$2; aE=$3; gS=$4; gE=$5;
        printf "%s/%s_%d_%d.fa\t%s_%d_%d_%d_%d\n", MONO, chr, gS, gE, chr, aS, aE, gS, gE;
      }
    ' "$bed" | while read -r mf header; do
      if [[ -f "$mf" ]]; then
        # rewrite FASTA header to include array start/end
        awk -v H="$header" '
          /^>/ {print ">" H; next} {print}
        ' "$mf"
      fi
    done > "$outfa"

    [ -s "$outfa" ] || rm -f "$outfa"
  done
done

# 4) Gather all per-bin FASTAs into one folder
BIN_FLAT_DIR="$OUTDIR/${PREFIX}_bin_monomers"
mkdir -p "$BIN_FLAT_DIR"
find "$OUTDIR" -type f -path "*_${PREFIX}_array_monomer/*.fa" -exec cp -f {} "$BIN_FLAT_DIR"/ \;
find "$OUTDIR" -type f -path "*_${PREFIX}_array_monomer/*.fasta" -exec cp -f {} "$BIN_FLAT_DIR"/ \;

echo "[✓] Done."
echo "  Combined bin table : $COMBINED_ALL"
echo "  Per-bin FASTAs     : <chr>_${PREFIX}_array_monomer/*.fa"
echo "  All bin FASTAs     : $BIN_FLAT_DIR/"
