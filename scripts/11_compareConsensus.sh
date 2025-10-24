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
# HiReNET: compareConsensus.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   BLAT self-comparisons of consensus HOR monomers across one or
#   many chromosomes and multiple identity thresholds.
#
# Behavior:
#   • If --chr has ONE item (e.g., "chr1"), BLATs per-threshold files:
#       <consensdir>/op_group/chr1_clust.cons.<thr>.fa  vs itself
#   • If --chr has MULTIPLE items (e.g., "chr1,chr2"), BLATs the
#       merged-per-threshold files:
#       <consensdir>/op_group/clust.cons.<thr>.fa       vs itself
#
# Input:
#   --chr CSV          One or more chromosomes (e.g., "chr1" or "chr1,chr2")
#   --consensdir DIR   Path to all_recluster_consensus_monomer/ or its
#                      op_group subfolder (auto-detected)
#   --outdir DIR       Output directory (created if missing)
#
# Dependencies:
#   - blat
#   - awk, sed, grep, xargs
#
# Output (under --outdir):
#   - blat_output_consensus/*.blat
#   - blat_output_consensus_sub/*.blat.sub
#   - merge_blat_consensus.txt
#
# Example:
#   HiReNET compareConsensus \
#       --chr "chr1,chr2" \
#       --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer \
#       --outdir AthCEN178_compare_consensusHOR_all
#
#   HiReNET compareConsensus \
#       --chr "chr1" \
#       --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer/op_group \
#       --outdir AthCEN178_compare_consensusHOR_chr1
#
# Notes:
#   - Thresholds default to 0.90–0.99 but can be overridden with --thr.
#   - If --consensdir points to the parent, the script auto-switches
#     to the 'op_group' subfolder when present.
# ==============================================================

die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing dependency in PATH: $1"; }

# defaults
CHR_CSV=""
CONSENSDIR=""
OUTDIR=""
THR_CSV="0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99"
MAXGAP=10
MINSCORE=0
REPMATCH=2147483647

usage() {
  cat <<'EOF'
============================================================
                  HiReNET: compareConsensus
============================================================

Description:
  BLAT self-comparisons of consensus HOR monomers. For a single
  chromosome, runs per-threshold self-comparisons. For multiple
  chromosomes, runs self-comparisons on the merged-per-threshold
  consensus files.

------------------------------------------------------------
Usage:
  HiReNET compareConsensus \
      --chr <csv> \
      --consensdir <path/to/all_recluster_consensus_monomer[ /op_group ]> \
      --outdir <path/to/output_dir> [options]

Required arguments:
  --chr <csv>            Chromosome list (e.g., "chr1" or "chr1,chr2")
  --consensdir <dir>     Parent consensus dir or its op_group subfolder
  --outdir <dir>         Output directory (created if missing)

Optional arguments:
  --thr <csv>            Identity thresholds (default: 0.90..0.99)
  --maxGap <int>         BLAT maxGap (default: 10)
  --minScore <int>       BLAT minScore (default: 0)
  --repMatch <int>       BLAT repMatch (default: 2147483647)
  -h, --help             Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • blat_output_consensus/*.blat
  • blat_output_consensus_sub/*.blat.sub
  • merge_blat_consensus.txt

Examples:
  HiReNET compareConsensus \
      --chr "chr1,chr2" \
      --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer \
      --outdir AthCEN178_compare_consensusHOR_all

  HiReNET compareConsensus \
      --chr "chr1" \
      --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer/op_group \
      --outdir AthCEN178_compare_consensusHOR_chr1

============================================================
EOF
}

# Backward compatibility (if older code calls print_help)
print_help() { usage; }

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

# deps
need blat
need awk

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
  # Compact table: match, qNameLen, tNameLen, qName, tName
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
