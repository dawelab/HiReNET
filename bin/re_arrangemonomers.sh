#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# re_arrangemonomers.sh 
#
# From a merged-bin classification table, collect monomers that fall
# inside each (chr, array, bin) and build per-bin FASTAs.
#
# Assumes --monomer-dir is a flat directory containing files named:
#   <chr>_<array_start>_<array_end>_<monomer_start>_<monomer_end>.fa
# e.g. chr1_14389697_18147662_100016_100193.fa
#
# Requires: awk, sed, grep, cut, sort, find, xargs, basename
# ============================================================

die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

print_help() {
  cat <<'EOF'
Usage:
  re_arrangemonomers.sh \
    --bins BINS_TSV \
    --class CLASS_NAME \
    --prefix PREFIX \
    --monomer-dir DIR_WITH_PER_MONOMER_FASTAS \
    --outdir OUTDIR \
    [--chr "chr1,chr2,..."] \
    [--chr-file chr_list.txt]

Notes:
  - BINS_TSV must have header with columns: chr, array, start, end
    and a class column among: lda_pred_smooth | lda_pred | class
  - monomer-dir must contain one FASTA per monomer, named:
      chr_arrayStart_arrayEnd_monomerStart_monomerEnd.fa
EOF
}

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

need awk; need sed; need grep; need cut; need sort; need find; need xargs; need basename

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

# ------------ filter bins to selected class ------------
BINS_BED="$OUTDIR/mergebin_bins_${CLASS}.bed"
awk -v FS='\t' -v OFS='\t' -v CLASS="$CLASS" -v CC="$CLASS_COL" '
  NR==1{
    for(i=1;i<=NF;i++) h[$i]=i
    next
  }
  $h[CC]==CLASS {
    # array is "start-end"
    n = split($h["array"], ab, /-/)
    if (n==2) print $h["chr"], ab[1]+0, ab[2]+0, $h["start"]+0, $h["end"]+0
  }
' "$BINS_FILE" > "$BINS_BED"
[[ -s "$BINS_BED" ]] || die "No bins matched class=$CLASS in $BINS_FILE"

# ------------ build per-chr monomer lists from flat dir ------------
echo "[*] Indexing per-chromosome monomers from filenames …"
for chr in "${CHRS[@]}"; do
  list_out="$OUTDIR/chr_list/${chr}_monomer_list.txt"
  # from filenames like: chr1_14389697_18147662_100016_100193.fa
  find "$MONOMER_DIR" -maxdepth 1 -type f -name "${chr}_*.fa" -printf '%f\n' \
    | sed 's/\.fa$//' \
    | awk -F'_' 'NF>=5 {print $1"\t"$2+0"\t"$3+0"\t"$4+0"\t"$5+0}' \
    | sort -k2,2n -k3,3n -k4,4n > "$list_out"
  [[ -s "$list_out" ]] || echo "[!] Empty monomer list for $chr at $list_out" >&2
done

# ------------ build per-bin lists + FASTAs ------------
echo "[*] Building per-bin lists and merged FASTAs …"
while read -r chr array_start array_end bin_start bin_end; do
  monolist="$OUTDIR/chr_list/${chr}_monomer_list.txt"
  listout="$OUTDIR/mergebin_list/${chr}_${array_start}_${array_end}_${bin_start}_${bin_end}.list.txt"
  [[ -s "$monolist" ]] || { echo "[!] Missing $monolist — skipping bin" >&2; continue; }

  # Pick monomers whose monomer_start falls strictly inside the bin
  awk -v C="$chr" -v A="$array_start" -v B="$array_end" -v S="$bin_start" -v E="$bin_end" '
    BEGIN{OFS="_"}
    $1==C && $2==A && $3==B && ($4+0) > S && ($4+0) < E {
      print $1,$2,$3,$4,$5".fa"
    }
  ' "$monolist" > "$listout"

  merged="$OUTDIR/re_arrange_monomers/${chr}_${array_start}_${array_end}_${bin_start}_${bin_end}.fa"
  : > "$merged"

  if [[ -s "$listout" ]]; then
    # concatenate files directly from flat monomer dir
    while read -r fn; do
      src="$MONOMER_DIR/$fn"
      if [[ -s "$src" ]]; then
        cat "$src" >> "$merged"
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
  base="$(basename "$fafile" .fa)"
  out="$BEDDIR/${base}.bed"
  bin_start="$(echo "$base" | awk -F'_' '{print $(NF-1)}')"
  bin_end="$(echo   "$base" | awk -F'_' '{print $(NF)}')"
  # Each record header is the monomer ID; reproduce 7 columns:
  # chr array_start array_end monomer_start monomer_end bin_start bin_end
  grep '^>' "$fafile" \
    | sed 's/^>//' \
    | awk -v OFS='\t' -v bs="$bin_start" -v be="$bin_end" -F'_' '
        NF>=5 {print $1, $2+0, $3+0, $4+0, $5+0, bs+0, be+0}
      ' > "$out"
done

COMBINED="$OUTDIR/${PREFIX}_monomer_bed_inbin.txt"
cat "$BEDDIR"/*.bed > "$COMBINED" 2>/dev/null || :

echo "[✓] Done."
echo "  Per-chr monomer lists : $OUTDIR/chr_list/"
echo "  Per-bin monomer lists : $OUTDIR/mergebin_list/"
echo "  Merged bin FASTAs     : $OUTDIR/re_arrange_monomers/"
echo "  Per-bin BEDs          : $OUTDIR/re_arrange_monomers_bed/"
echo "  Combined BED          : $COMBINED"