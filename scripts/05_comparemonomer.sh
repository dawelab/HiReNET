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
# HiReNET: comparemonomer.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Compare per-bin monomer FASTAs using BLAT to evaluate monomer
#   sequence similarity within and between bins. Performs:
#     1) All-to-all self-BLAT per bin (always)
#     2) Bin-vs-consensus BLAT (optional, if --consensus provided)
#
# Dependencies:
#   - blat
#   - awk
#   - sed
#   - xargs
#
# Input:
#   --bins-dir     Directory containing per-bin FASTAs (*.fa|*.fasta)
#   --consensus    Optional consensus FASTA for bin-vs-consensus BLAT
#
# Output (under --outdir):
#   - blat_output/*.blat
#   - blat_output_sub/*.blat.sub
#   - merge_blat.txt
#   - [optional] blat_con_output/*.con.blat
#   - [optional] blat_con_output_sub/*.con.blat.sub
#   - [optional] merge_blat_con.txt
#
# Example:
#   HiReNET comparemonomer \
#       --bins-dir AthCEN178_arrangemonomer_10kb/AthCEN178_bin_monomers \
#       --outdir   AthCEN178_comparemonomers \
#       --consensus AthCEN178_consensus.fasta
#
# Notes:
#   - Runs BLAT self-comparisons for every bin FASTA.
#   - Optional consensus comparison adds cross-bin validation.
#   - Outputs are BLAT tables for downstream classification.
# ==============================================================

# ---- defaults ----
BINS_DIR=""
OUTDIR=""
CONSENSUS=""     # optional; if empty, skip consensus step
MAXGAP=10
MINSCORE=0
REPMATCH=2147483647

die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

usage() {
  cat <<'EOF'
============================================================
                     HiReNET: comparemonomer
============================================================

Description:
  Perform BLAT-based comparisons of per-bin monomer FASTAs.
  Always runs self-BLAT (each bin compared to itself) and,
  optionally, compares bins to consensus sequences if provided.

------------------------------------------------------------
Usage:
  HiReNET comparemonomer \
      --bins-dir <path/to/bin_fastas> \
      --outdir   <path/to/output_dir> [options]

Required arguments:
  --bins-dir <dir>        Directory containing per-bin FASTAs (*.fa|*.fasta)
  --outdir <dir>          Output directory (created if missing)

Optional arguments:
  --consensus <file>      FASTA of consensus sequences; if provided, also run
                          bin-vs-consensus BLAT and generate:
                            outdir/blat_con_output/*.con.blat
                            outdir/blat_con_output_sub/*.con.blat.sub
                            outdir/merge_blat_con.txt
  --maxGap <int>          BLAT maxGap (default: 10)
  --minScore <int>        BLAT minScore (default: 0)
  --repMatch <int>        BLAT repMatch (default: 2147483647)
  -h, --help              Show this help message and exit

------------------------------------------------------------
Outputs (in <outdir>):
  • blat_output/*.blat
  • blat_output_sub/*.blat.sub
  • merge_blat.txt
  • [optional] blat_con_output/*.con.blat
  • [optional] blat_con_output_sub/*.con.blat.sub
  • [optional] merge_blat_con.txt

Example:
  HiReNET comparemonomer \
      --bins-dir AthCEN178_arrangemonomer_10kb/AthCEN178_bin_monomers \
      --outdir   AthCEN178_comparemonomers \
      --consensus AthCEN178_consensus.fasta

============================================================
EOF
}

# Backward compatibility (if older calls use print_help)
print_help() { usage; }

[[ $# -eq 0 ]] && { print_help; exit 1; }

# ---- parse args ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bins-dir)    BINS_DIR="${2:-}"; shift 2;;
    --outdir)      OUTDIR="${2:-}"; shift 2;;
    --consensus)   CONSENSUS="${2:-}"; shift 2;;
    --maxGap)      MAXGAP="${2:-}"; shift 2;;
    --minScore)    MINSCORE="${2:-}"; shift 2;;
    --repMatch)    REPMATCH="${2:-}"; shift 2;;
    -h|--help)     print_help; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -z "$BINS_DIR" ]] && die "Provide --bins-dir"
[[ -z "$OUTDIR"   ]] && die "Provide --outdir"

BINS_DIR="$(readlink -f "$BINS_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"
mkdir -p "$OUTDIR"

# Resolve consensus only if provided
if [[ -n "$CONSENSUS" ]]; then
  CONSENSUS="$(readlink -f "$CONSENSUS")"
fi

# ---- deps (PATH-based) ----
need blat
need awk
need sed
need xargs

echo "[*] bins-dir:   $BINS_DIR"
echo "[*] outdir:     $OUTDIR"
if [[ -n "$CONSENSUS" ]]; then
  echo "[*] consensus:  $CONSENSUS"
else
  echo "[*] consensus:  (none) — will skip consensus BLAT"
fi
echo "[*] BLAT opts:  -maxGap=$MAXGAP -minScore=$MINSCORE -repMatch=$REPMATCH"

# ---- gather bin FASTAs ----
shopt -s nullglob
BIN_FASTAS=( "$BINS_DIR"/*.fa "$BINS_DIR"/*.fasta )
((${#BIN_FASTAS[@]})) || die "No per-bin FASTAs (*.fa|*.fasta) found in $BINS_DIR"

# Create required dirs
mkdir -p "$OUTDIR/blat_output" "$OUTDIR/blat_output_sub"
if [[ -n "$CONSENSUS" ]]; then
  mkdir -p "$OUTDIR/blat_con_output" "$OUTDIR/blat_con_output_sub"
fi

# =========================================================
# 1) ALL-TO-ALL (self-BLAT per bin)  [always]
# =========================================================
echo "[*] Self-BLAT (all-to-all) on ${#BIN_FASTAS[@]} bins…"
for bin in "${BIN_FASTAS[@]}"; do
  name="$(basename "$bin")"
  stem="${name%.*}"                    # drop extension
  out_blat="$OUTDIR/blat_output/${stem}.blat"

  echo "  [-] BLAT self: $name"
  blat "$bin" "$bin" \
       -t=dna -q=dna \
       -maxGap="$MAXGAP" -minScore="$MINSCORE" -repMatch="$REPMATCH" \
       "$out_blat"

  # Compact table
  awk '$1 ~ /^[0-9]+$/' "$out_blat" \
    | awk 'NR>1' \
    | awk '{print $1"\t"$10"\t"$11"\t"$14"\t"$15}' \
    > "$OUTDIR/blat_output_sub/${stem}.blat.sub"
done

# =========================================================
# 2) ALL-TO-CONSENSUS (bin vs consensus)  [only if set]
# =========================================================
if [[ -n "$CONSENSUS" ]]; then
  echo "[*] BLAT to consensus on ${#BIN_FASTAS[@]} bins…"
  for bin in "${BIN_FASTAS[@]}"; do
    name="$(basename "$bin")"
    stem="${name%.*}"
    out_blat="$OUTDIR/blat_con_output/${stem}.con.blat"

    echo "  [-] BLAT vs consensus: $name"
    blat "$bin" "$CONSENSUS" \
         -t=dna -q=dna \
         -maxGap="$MAXGAP" -minScore="$MINSCORE" -repMatch="$REPMATCH" \
         "$out_blat"

    awk '$1 ~ /^[0-9]+$/' "$out_blat" \
      | awk 'NR>1' \
      | awk '{print $1"\t"$10"\t"$11"\t"$14"\t"$15}' \
      > "$OUTDIR/blat_con_output_sub/${stem}.con.blat.sub"
  done

  echo "[*] Merging consensus subtables…"
  MERGE="$OUTDIR/merge_blat_con.txt"
  : > "$MERGE"
  for sub in "$OUTDIR"/blat_con_output_sub/*.con.blat.sub; do
    [[ -e "$sub" ]] || break
    cat "$sub" >> "$MERGE"
  done
fi

echo "[✓] Done."
echo "Outputs:"
echo "  $OUTDIR/blat_output/*.blat"
echo "  $OUTDIR/blat_output_sub/*.blat.sub"
if [[ -n "$CONSENSUS" ]]; then
  echo "  $OUTDIR/blat_con_output/*.con.blat"
  echo "  $OUTDIR/blat_con_output_sub/*.con.blat.sub"
  echo "  $OUTDIR/merge_blat_con.txt"
else
  echo "  (consensus outputs skipped; no --consensus provided)"
fi
