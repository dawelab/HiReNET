#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# arrangeHORmonomer.sh
# 1) Create --outdir
# 2) Copy GROUPDIR/group_HOR_monomers -> OUTDIR/group_HOR_monomers
# 3) Build OUTDIR/group_HOR_monomers_output by concatenating monomer FASTAs
#
# Usage:
#   ./arrangeHORmonomer.sh \
#     --groupdir network_mergebin \
#     --monomer-dir monomerout_partgenome/AthCEN178_monomers \
#     --outdir network_mergebin_consensus
#
# Requires: awk, xargs, rsync (or cp -a), grep
# ------------------------------------------------------------

GROUPDIR=""
OUTDIR=""
MONOMER_DIR=""

usage() {
  cat <<'USAGE'
arrangeHORmonomer.sh
Required:
  --groupdir        Path containing group_HOR_monomers/* (per-bin folders with *.list.txt)
  --monomer-dir    Flat folder of single-monomer FASTAs (names listed in *.list.txt)
  --outdir         Output base (will be created). Copies group_HOR_monomers here.

Example:
  ./arrangeHORmonomer.sh \
    --groupdir network_mergebin \
    --monomer-dir monomerout_partgenome/AthCEN178_monomers \
    --outdir network_mergebin_consensus
USAGE
  exit 1
}

[[ $# -eq 0 ]] && usage
while [[ $# -gt 0 ]]; do
  case "$1" in
    --groupdir)       GROUPDIR="${2:-}"; shift 2;;
    --groupdir=*)     GROUPDIR="${1#*=}"; shift;;
    --monomer-dir)   MONOMER_DIR="${2:-}"; shift 2;;
    --monomer-dir=*) MONOMER_DIR="${1#*=}"; shift;;
    --outdir)        OUTDIR="${2:-}"; shift 2;;
    --outdir=*)      OUTDIR="${1#*=}"; shift;;
    -h|--help)       usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -z "$GROUPDIR"     ]] && { echo "ERROR: --groupdir is required"; usage; }
[[ -z "$MONOMER_DIR" ]] && { echo "ERROR: --monomer-dir is required"; usage; }
[[ -z "$OUTDIR"      ]] && { echo "ERROR: --outdir is required"; usage; }

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing <$1>" >&2; exit 2; }; }
need awk; need xargs; need grep
if ! command -v rsync >/dev/null 2>&1; then
  echo "WARN: rsync not found; falling back to cp -a" >&2
fi

GROUPDIR="$(readlink -f "$GROUPDIR")"
MONOMER_DIR="$(readlink -f "$MONOMER_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"

mkdir -p "$OUTDIR"

src_lists="$GROUPDIR/group_HOR_monomers"
dst_lists="$OUTDIR/group_HOR_monomers"
[[ -d "$src_lists" ]] || { echo "ERROR: missing $src_lists"; exit 3; }

# Copy lists into OUTDIR
if command -v rsync >/dev/null 2>&1; then
  rsync -a --delete "$src_lists"/ "$dst_lists"/
else
  mkdir -p "$dst_lists"
  (shopt -s dotglob; cp -a "$src_lists"/* "$dst_lists"/)
fi

# Build group_HOR_monomers_output
out_lists="$OUTDIR/group_HOR_monomers_output"
mkdir -p "$out_lists"

shopt -s nullglob
for d in "$dst_lists"/*; do
  [[ -d "$d" ]] || continue
  dirbase="$(basename "$d")"
  outbin="$out_lists/$dirbase"
  mkdir -p "$outbin"

  for file in "$d"/*.list.txt; do
    [[ -e "$file" ]] || continue
    name="$(basename "$file" .list.txt)"
    outfile="$outbin/$name.fa"
    awk -v m="$MONOMER_DIR" '{print m "/" $0}' "$file" | xargs -r cat -- > "$outfile"
  done
done

echo "[✓] arrangeHORmonomer: done"
echo "  Copied lists      : $dst_lists/"
echo "  Re-clustered FASTAs -> $out_lists/"