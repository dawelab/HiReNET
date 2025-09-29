#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# compareConsensus.sh
# BLAT consensus monomers in two modes:
#   --chr "chr1"            ➜ BLAT   <consensdir>/op_group/chr1_clust.cons.91.fa ... .99.fa (self)
#   --chr "chr1,chr2,..."   ➜ BLAT   <consensdir>/op_group/clust.cons.91.fa   ... .99.fa (self)
#
# You may pass --consensdir either as the parent dir (.../all_recluster_consensus_monomer)
# or directly as the op_group dir (.../all_recluster_consensus_monomer/op_group).
# The script will auto-switch to the 'op_group' subfolder if present.
#
# Usage examples:
#   ./compareConsensus.sh --chr "chr1,chr2" \
#     --consensdir /scratch/.../all_recluster_consensus_monomer \
#     --outdir compare_consensusHOR_all
#
#   ./compareConsensus.sh --chr "chr1" \
#     --consensdir /scratch/.../all_recluster_consensus_monomer/op_group \
#     --outdir compare_consensusHOR_all
# =========================================================

die(){ echo "ERROR: $*" >&2; exit 1; }

# defaults
CHR_CSV=""
CONSENSDIR=""
OUTDIR=""
THR_CSV="0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99"
MAXGAP=10
MINSCORE=0
REPMATCH=2147483647

print_help(){
  cat <<'EOF'
Usage:
  compareConsensus.sh --chr CSV --consensdir DIR --outdir DIR [--thr CSV] [--maxGap N] [--minScore N] [--repMatch N]

Notes:
  - If --chr has ONE item (e.g. "chr1"), we BLAT per-threshold files:
      <consensdir>/op_group/<chr>_clust.cons.<thr>.fa  vs itself
  - If --chr has MULTIPLE items (e.g. "chr1,chr2"), we BLAT merged-per-threshold files:
      <consensdir>/op_group/clust.cons.<thr>.fa        vs itself
  - --consensdir can be the parent dir or the op_group dir; op_group is auto-detected.
EOF
}

[[ $# -eq 0 ]] && { print_help; exit 1; }

# parse
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chr)        CHR_CSV="${2:-}"; shift 2;;
    --consensdir) CONSENSDIR="${2:-}"; shift 2;;
    --outdir)     OUTDIR="${2:-}"; shift 2;;
    --thr)        THR_CSV="${2:-}"; shift 2;;
    --maxGap)     MAXGAP="${2:-}"; shift 2;;
    --minScore)   MINSCORE="${2:-}"; shift 2;;
    --repMatch)   REPMATCH="${2:-}"; shift 2;;
    -h|--help)    print_help; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

[[ -z "$CHR_CSV"    ]] && die "Provide --chr"
[[ -z "$CONSENSDIR" ]] && die "Provide --consensdir"
[[ -z "$OUTDIR"     ]] && die "Provide --outdir"

CONSENSDIR="$(readlink -f "$CONSENSDIR")"
# auto-detect op_group
if [[ -d "$CONSENSDIR/op_group" ]]; then
  CONSENSDIR="$CONSENSDIR/op_group"
fi
[[ -d "$CONSENSDIR" ]] || die "Consensus dir not found: $CONSENSDIR"

OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"/{blat,blat_sub}

IFS=',' read -r -a CHR_ARR <<< "$CHR_CSV"
IFS=',' read -r -a THR_ARR <<< "$THR_CSV"

# BLAT availability
if command -v ml >/dev/null 2>&1; then ml BLAT/3.7-GCC-12.3.0 || true; fi
command -v blat >/dev/null 2>&1 || die "blat not found in PATH"

echo "[*] consensdir : $CONSENSDIR"
echo "[*] outdir     : $OUTDIR"
echo "[*] chr list   : ${CHR_ARR[*]}"
echo "[*] thresholds : ${THR_ARR[*]}"
echo "[*] BLAT opts  : -maxGap=$MAXGAP -minScore=$MINSCORE -repMatch=$REPMATCH"

BLAT_OUT="$OUTDIR/blat"
BLAT_SUB="$OUTDIR/blat_sub"

run_blat_self(){
  local src="$1" base="$2"
  local out_blat="$BLAT_OUT/${base}.blat"
  local out_sub="$BLAT_SUB/${base}.blat.sub"
  [[ -s "$src" ]] || { echo "[!] Missing/empty: $src"; return; }
  echo "  [-] BLAT self: $(basename "$src")"
  blat "$src" "$src" -t=dna -q=dna \
       -maxGap="$MAXGAP" -minScore="$MINSCORE" -repMatch="$REPMATCH" \
       "$out_blat"
  awk '$1 ~ /^[0-9]+$/ && NR>1 {print $1"\t"$10"\t"$11"\t"$14"\t"$15}' \
    "$out_blat" > "$out_sub"
}

if [[ ${#CHR_ARR[@]} -le 1 ]]; then
  # single-chromosome mode
  CHR="${CHR_ARR[0]}"
  echo "[*] Mode: single chromosome ($CHR)"
  for th in "${THR_ARR[@]}"; do
    th_out="${th/0./}"  # 0.91 -> 91
    src="$CONSENSDIR/${CHR}_clust.cons.${th_out}.fa"
    base="${CHR}.cons.${th_out}"
    run_blat_self "$src" "$base"
  done
else
  # merged-per-threshold mode
  echo "[*] Mode: merged per-threshold (>=2 chromosomes)"
  for th in "${THR_ARR[@]}"; do
    th_out="${th/0./}"
    src="$CONSENSDIR/clust.cons.${th_out}.fa"
    base="merged.cons.${th_out}"
    run_blat_self "$src" "$base"
  done
fi

echo "[✓] Done."
echo "  BLAT outputs : $BLAT_OUT/"
echo "  Subtables    : $BLAT_SUB/"