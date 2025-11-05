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
# HiReNET: arrangeHORmonomer.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Assemble and organize HOR monomer groups for consensus building.  
#   This step collects monomer FASTAs for each HOR cluster group 
#   (generated from networkHOR) and prepares merged FASTA files 
#   for downstream consensus generation.
#
# Workflow:
#   1) Create output directory (--outdir)
#   2) Copy group_HOR_monomers/ from --groupdir
#   3) Concatenate monomer FASTAs into:
#        --outdir/group_HOR_monomers_output/
#
# Input:
#   --groupdir <dir>     Directory containing group_HOR_monomers/* 
#                        (each subfolder includes *.list.txt with monomer names)
#   --monomer-dir <dir>  Directory containing individual monomer FASTAs
#                        named like chr_arrayStart_arrayEnd_monoStart_monoEnd.fa
#   --outdir <dir>       Output directory (created if missing)
#
# Dependencies:
#   - awk, grep, xargs, cp, sed, rsync (optional)
#
# Output (under --outdir):
#   - group_HOR_monomers/
#   - group_HOR_monomers_output/
#       ├── groupX_combined.fa
#       ├── groupY_combined.fa
#
# Example:
#   HiReNET arrangeHORmonomer \
#       --groupdir AthCEN178_network_HOR_mergebin \
#       --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
#       --outdir AthCEN178_network_mergebin_consensus
#
# Notes:
#   - Each *.list.txt file inside group_HOR_monomers specifies 
#     which monomer FASTAs to concatenate.
#   - The resulting combined FASTAs serve as input for 
#     HiReNET consensusHORmonomer (step 10).
# ==============================================================

GROUPDIR=""
OUTDIR=""
MONOMER_DIR=""

usage() {
  cat <<'EOF'
============================================================
                 HiReNET: arrangeHORmonomer
============================================================

Description:
  Gather and combine HOR monomer FASTAs for each HOR cluster.  
  This step organizes outputs from networkHOR into grouped 
  monomer files, preparing them for consensus sequence generation.

------------------------------------------------------------
Usage:
  HiReNET arrangeHORmonomer \
      --groupdir <path/to/networkHOR_out> \
      --monomer-dir <path/to/monomer_fastas> \
      --outdir <path/to/output_dir>

Required arguments:
  --groupdir <dir>       Directory containing group_HOR_monomers/*
  --monomer-dir <dir>    Directory of individual monomer FASTAs
  --outdir <dir>         Output directory (created if missing)

------------------------------------------------------------
Outputs (in <outdir>):
  • group_HOR_monomers/  
  • group_HOR_monomers_output/
      ├── groupX_combined.fa  
      ├── groupY_combined.fa  

Example:
  HiReNET arrangeHORmonomer \
      --groupdir AthCEN178_network_HOR_mergebin \
      --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
      --outdir AthCEN178_network_mergebin_consensus

============================================================
EOF
}

# Backward compatibility (for older calls)
print_help() { usage; }

die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing <$1>" >&2; exit 2; }; }

[[ $# -eq 0 ]] && usage
while [[ $# -gt 0 ]]; do
  case "$1" in
    --groupdir)       GROUPDIR="${2:-}"; shift 2;;
    --groupdir=*)     GROUPDIR="${1#*=}"; shift;;
    --monomer-dir)    MONOMER_DIR="${2:-}"; shift 2;;
    --monomer-dir=*)  MONOMER_DIR="${1#*=}"; shift;;
    --outdir)         OUTDIR="${2:-}"; shift 2;;
    --outdir=*)       OUTDIR="${1#*=}"; shift;;
    -h|--help)        usage;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -z "$GROUPDIR"    ]] && { echo "ERROR: --groupdir is required"; usage; }
[[ -z "$MONOMER_DIR" ]] && { echo "ERROR: --monomer-dir is required"; usage; }
[[ -z "$OUTDIR"      ]] && { echo "ERROR: --outdir is required"; usage; }

need awk; need xargs; need grep; need sed
# rsync is optional (falls back to cp -a if absent)

GROUPDIR="$(readlink -f "$GROUPDIR")"
MONOMER_DIR="$(readlink -f "$MONOMER_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"

mkdir -p "$OUTDIR"

src_lists="$GROUPDIR/group_HOR_monomers"
dst_lists="$OUTDIR/group_HOR_monomers"
[[ -d "$src_lists" ]] || die "Missing $src_lists"

# Copy lists into OUTDIR
if command -v rsync >/dev/null 2>&1; then
  rsync -a --delete "$src_lists"/ "$dst_lists"/
else
  echo "WARN: rsync not found; using cp -a" >&2
  mkdir -p "$dst_lists"
  shopt -s dotglob
  cp -a "$src_lists"/* "$dst_lists"/ 2>/dev/null || true
  shopt -u dotglob
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
    : > "$outfile"

    # Each line in $file is chr_aS_aE_gS_gE.fa (or with a path prefix)
    while IFS= read -r line || [[ -n "$line" ]]; do
      # strip trailing CR and directory/path
      clean="$(echo "$line" | sed 's/\r$//' | sed 's#.*/##')"
      # bake out ".fa" or ".fasta"
      stem="$(echo "$clean" | sed 's/\.\(fa\|fasta\)$//')"

      # Expect fields: chr_aS_aE_gS_gE (preferred) OR chr_gS_gE (fallback)
      IFS='_' read -r f1 f2 f3 f4 f5 <<< "$stem"

      if [[ -n "${f5:-}" ]]; then
        chr="$f1"; aS="$f2"; aE="$f3"; gS="$f4"; gE="$f5"
      else
        # fallback: only chr_gS_gE in list; use aS/aE from filename later
        chr="$f1"; aS=""; aE=""; gS="$f2"; gE="$f3"
      fi

      src1="$MONOMER_DIR/${chr}_${gS}_${gE}.fa"
      src2="$MONOMER_DIR/${chr}_${gS}_${gE}.fasta"
      src=""
      [[ -s "$src1" ]] && src="$src1"
      [[ -z "$src" && -s "$src2" ]] && src="$src2"

      if [[ -n "$src" ]]; then
        header="${chr}_${aS:-NA}_${aE:-NA}_${gS}_${gE}"
        awk -v H="$header" 'BEGIN{printed=0}
          /^>/ {print ">" H; next}
               {print}
        ' "$src" >> "$outfile"
      else
        echo "[!] Missing monomer FASTA for ${chr}_${gS}_${gE} in $MONOMER_DIR" >&2
      fi
    done < "$file"

    # drop empty outputs
    [[ -s "$outfile" ]] || rm -f "$outfile"
  done
done
shopt -u nullglob

echo "[✓] arrangeHORmonomer: done"
echo "  Copied lists          : $dst_lists/"
echo "  Re-clustered FASTAs   : $out_lists/"
