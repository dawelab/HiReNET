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
# HiReNET: rearrangemonomers.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Reconstruct and group monomer FASTAs according to HOR class 
#   predictions across merged genomic bins. For each (chr, array, 
#   bin, class), monomers are gathered into per-bin FASTAs to 
#   enable downstream HOR pattern analysis and consensus building.
#
# Input:
#   --bins <file>          Classification table with columns:
#                            chr, array, start, end
#                          where:
#                            • array = "arrayStart-arrayEnd" (array coordinates)
#                            • start/end = genomic bin boundaries
#
#   --class <string>       Class name or label to extract (e.g., HOR1, HOR2)
#   --prefix <string>      Prefix for labeling all output files
#   --monomer-dir <dir>    Directory with per-monomer FASTAs named as:
#                            chr_arrayStart_arrayEnd_monoStart_monoEnd.fa
#                          (monomerStart/End are array-relative coordinates)
#   --outdir <dir>         Output directory for grouped FASTAs
#
# Dependencies:
#   - awk, sed, grep, cut, sort, find, xargs, basename
#
# Output (under --outdir):
#   - ${PREFIX}_${CLASS}_bins/
#       ├── chrX_arrayStart_arrayEnd_binStart_binEnd.fa
#       ├── combined_${PREFIX}_${CLASS}.fa
#   - Summary table of grouped monomers per bin
#
# Example:
#   HiReNET rearrangemonomers \
#       --bins AthCEN178_classpred_out/AthCEN178_bin_class_smooth.txt \
#       --class HOR1 \
#       --prefix AthCEN178 \
#       --monomer-dir AthCEN178_arrangemonomer_10kb/AthCEN178_bin_monomers \
#       --outdir AthCEN178_rearranged_HOR1
#
# Notes:
#   - The bins file must correspond to smoothed outputs from 
#     HiReNET classprediction (step 06).
#   - Monomer FASTA names must exactly match HiReNET’s naming scheme.
# ==============================================================

die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

usage() {
  cat <<'EOF'
============================================================
                   HiReNET: rearrangemonomers
============================================================

Description:
  Reorganize monomer FASTAs into per-bin groups based on 
  classification output tables. Each bin and HOR class will 
  have its own FASTA for downstream HOR pattern analysis.

------------------------------------------------------------
Usage:
  HiReNET rearrangemonomers \
      --bins <path/to/bin_class_table.txt> \
      --class <class_name> \
      --prefix <prefix> \
      --monomer-dir <path/to/monomer_fastas> \
      --outdir <path/to/output_dir> [options]

Required arguments:
  --bins <file>          Classification table (e.g., *_bin_class_smooth.txt)
  --class <string>       Class label to extract (e.g., HOR1, HOR2)
  --prefix <string>      Prefix for naming output files
  --monomer-dir <dir>    Directory with per-monomer FASTAs
  --outdir <dir>         Output directory (created if missing)

Optional arguments:
  --chr "chr1,chr2,..."  Comma-separated list of chromosomes
  --chr-file <file>      File containing chromosome names (one per line)
  -h, --help             Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • ${PREFIX}_${CLASS}_bins/
      ├── chrX_arrayStart_arrayEnd_binStart_binEnd.fa
      ├── combined_${PREFIX}_${CLASS}.fa
  • Summary table of monomers grouped per bin

Example:
  HiReNET rearrangemonomers \
      --bins AthCEN178_classpred_out/AthCEN178_bin_class_smooth.txt \
      --class HOR1 \
      --prefix AthCEN178 \
      --monomer-dir AthCEN178_arrangemonomer_10kb/AthCEN178_bin_monomers \
      --outdir AthCEN178_rearranged_HOR1

============================================================
EOF
}

# Backward compatibility (if older scripts still call print_help)
print_help() { usage; }

[[ $# -eq 0 ]] && { print_help; exit 1; }

# ------------ args -----------
BINS_FILE=""
CLASS=""
MONOMER_DIR=""
PREFIX=""
OUTDIR=""
CHR_CSV=""
CHR_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bins)         BINS_FILE="$2"; shift 2;;
    --class)        CLASS="$2"; shift 2;;
    --monomer-dir)  MONOMER_DIR="$2"; shift 2;;
    --prefix)       PREFIX="$2"; shift 2;;
    --outdir)       OUTDIR="$2"; shift 2;;
    --chr)          CHR_CSV="$2"; shift 2;;
    --chr-file)     CHR_FILE="$2"; shift 2;;
    -h|--help)      print_help; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "$BINS_FILE"   ]] && die "Provide --bins"
[[ -z "$CLASS"       ]] && die "Provide --class"
[[ -z "$MONOMER_DIR" ]] && die "Provide --monomer-dir"
[[ -z "$PREFIX"      ]] && die "Provide --prefix"
[[ -z "$OUTDIR"      ]] && die "Provide --outdir"

# ------------ deps ------------
need awk; need sed; need grep; need cut; need sort; need find; need xargs; need basename

# ------------ paths -----------
BINS_FILE="$(readlink -f "$BINS_FILE")"
MONOMER_DIR="$(readlink -f "$MONOMER_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"

mkdir -p "$OUTDIR"/{chr_list,mergebin_list,re_arrange_monomers,re_arrange_monomers_bed}

echo "[*] bins-file     : $BINS_FILE"
echo "[*] class filter  : $CLASS"
echo "[*] monomer-dir   : $MONOMER_DIR"
echo "[*] prefix        : $PREFIX"
echo "[*] outdir        : $OUTDIR"

# ------------ choose chromosomes ------------
declare -a CHRS
if [[ -n "$CHR_FILE" ]]; then
  mapfile -t CHRS < <(grep -v '^[[:space:]]*$' "$CHR_FILE")
elif [[ -n "$CHR_CSV" ]]; then
  IFS=',' read -r -a CHRS <<< "$CHR_CSV"
else
  # infer from monomer filenames present
  mapfile -t CHRS < <(find "$MONOMER_DIR" -maxdepth 1 -type f -name '*.fa' \
                        -printf '%f\n' \
                        | awk -F'_' '{print $1}' \
                        | sort -uV)
fi
((${#CHRS[@]})) || die "No chromosomes resolved (use --chr/--chr-file, or ensure monomer FASTAs exist)"

echo "[*] chromosomes   : ${CHRS[*]}"

# ------------ validate bins file -------------
header="$(head -n1 "$BINS_FILE" | tr -d '\r')"
for col in chr array start end; do
  echo "$header" | grep -wq "$col" || die "Bins file must have header column: $col"
done
CLASS_COL=""
for c in lda_pred_smooth lda_pred class; do
  if echo "$header" | grep -wq "$c"; then CLASS_COL="$c"; break; fi
done
[[ -z "$CLASS_COL" ]] && die "Bins file must include one of: lda_pred_smooth | lda_pred | class"
echo "[*] class column  : $CLASS_COL"

# ------------ filter to class; write chr aS aE bStart(genomic) bEnd(genomic) ------------
BINS_BED="$OUTDIR/mergebin_bins_${CLASS}.bed"
awk -v FS='\t' -v OFS='\t' -v CLASS="$CLASS" -v CC="$CLASS_COL" '
  NR==1{
    for(i=1;i<=NF;i++) h[$i]=i
    next
  }
  $h[CC]==CLASS {
    n = split($h["array"], ab, /-/)
    if (n==2) {
      chr = $h["chr"]
      aS  = ab[1]+0
      aE  = ab[2]+0
      bS  = $h["start"]+0
      bE  = $h["end"]+0
      print chr, aS, aE, bS, bE
    }
  }
' "$BINS_FILE" > "$BINS_BED"
[[ -s "$BINS_BED" ]] || die "No bins matched class=$CLASS in $BINS_FILE"

# ------------ build per-chr monomer lists from flat dir ------------
echo "[*] Indexing per-chromosome monomers from filenames …"
for chr in "${CHRS[@]}"; do
  list_out="$OUTDIR/chr_list/${chr}_monomer_list.txt"
  # Columns emitted: chr  gS  gE  filename
  find "$MONOMER_DIR" -maxdepth 1 -type f -name "${chr}_*.fa" -printf '%f\n' \
    | sed 's/\.fa$//' \
    | awk -F'_' 'NF==3 {printf "%s\t%d\t%d\t%s.fa\n", $1, $2+0, $3+0, $0}' \
    | sort -k2,2n -k3,3n > "$list_out"
  [[ -s "$list_out" ]] || echo "[!] Empty monomer list for $chr at $list_out" >&2
done

# ------------ build per-bin lists + FASTAs ------------
echo "[*] Building per-bin lists and merged FASTAs …"
while read -r chr array_start array_end bin_start bin_end; do
  monolist="$OUTDIR/chr_list/${chr}_monomer_list.txt"
  listout="$OUTDIR/mergebin_list/${chr}_${array_start}_${array_end}_${bin_start}_${bin_end}.list.txt"
  [[ -s "$monolist" ]] || { echo "[!] Missing $monolist — skipping bin" >&2; continue; }

  # pick monomers whose genomic start falls in [bin_start, bin_end)
  awk -v C="$chr" -v S="$bin_start" -v E="$bin_end" \
      '$1==C && ($2>=S) && ($2<E) {print $4}' "$monolist" > "$listout"

  merged="$OUTDIR/re_arrange_monomers/${chr}_${array_start}_${array_end}_${bin_start}_${bin_end}.fa"
  : > "$merged"

  if [[ -s "$listout" ]]; then
    while read -r fn; do
      src="$MONOMER_DIR/$fn"
      if [[ -s "$src" ]]; then
        # Rewrite header: >chr_gS_gE  →  >chr_aS_aE_gS_gE
        awk -v aS="$array_start" -v aE="$array_end" '
          /^>/ {
            sub(/^>/,"");
            n=$0; split(n, p, "_");  # p[1]=chr, p[2]=gS, p[3]=gE
            printf(">%s_%d_%d_%s_%s\n", p[1], aS, aE, p[2], p[3]);
            next
          }
          { print }
        ' "$src" >> "$merged"
      else
        echo "[!] Missing monomer FASTA: $src" >&2
      fi
    done < "$listout"
  fi
done < "$BINS_BED"

# ------------ per-bin BEDs ------------
echo "[*] Emitting per-bin BEDs …"
BEDDIR="$OUTDIR/re_arrange_monomers_bed"
mkdir -p "$BEDDIR"

for fafile in "$OUTDIR"/re_arrange_monomers/*.fa; do
  [[ -e "$fafile" ]] || continue

  base="$(basename "$fafile" .fa)"  # chr_aS_aE_bS_bE
  out="$BEDDIR/${base}.bed"

  # parse aS, aE, bS, bE from the filename
  aS=$(echo "$base" | awk -F'_' '{print $2}')
  aE=$(echo "$base" | awk -F'_' '{print $3}')
  bS=$(echo "$base" | awk -F'_' '{print $4}')
  bE=$(echo "$base" | awk -F'_' '{print $5}')

  # headers are now >chr_aS_aE_gS_gE (genomic), but support old >chr_gS_gE too
  grep '^>' "$fafile" \
    | sed 's/^>//' \
    | awk -v OFS="\t" -v aSfile="$aS" -v aEfile="$aE" -v bS="$bS" -v bE="$bE" -F'_' '
        NF==5 {
          chr=$1; aS=$2+0; aE=$3+0; gS=$4+0; gE=$5+0;
          print chr, aS, aE, gS, gE, bS+0, bE+0; next
        }
        NF==3 {
          chr=$1; gS=$2+0; gE=$3+0;
          print chr, aSfile+0, aEfile+0, gS, gE, bS+0, bE+0; next
        }
      ' > "$out"
done

COMBINED="$OUTDIR/${PREFIX}_monomer_bed_inbin.txt"
cat "$BEDDIR"/*.bed > "$COMBINED" 2>/dev/null || :
echo "[*] Wrote: $COMBINED"
