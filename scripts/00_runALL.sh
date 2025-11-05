#!/usr/bin/env bash
set -euo pipefail

# --- find HiReNET root, keep current working dir ---
if [[ -n "${HIRENET_ROOT:-}" && -d "${HIRENET_ROOT}" ]]; then
  ROOT="${HIRENET_ROOT}"
else
  ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
export HIRENET_ROOT="${ROOT}"

die(){ echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat <<'EOF'
============================================================
                       HiReNET: runALL
============================================================

Description:
  End-to-end runner: build HMM, find arrays, detect monomers,
  arrange into bins, compare, classify HORs, rearrange, network,
  consensus per chromosome, and shared-HOR plots.

------------------------------------------------------------
Usage:
  HiReNET runALL \
    --chr "chr1,chr2,..." \
    --prefix NAME \
    --consensus <consensus.fasta> \
    --seq <genome.fasta> \
    --variant <consensus_variant.fasta> \
    [--bin 10000] [--threads 10] [--plotv V2] [--min-monomer-len 120]

Required arguments:
  --chr         Comma-separated chromosomes (e.g., "chr1,chr2,chr3,chr4,chr5")
  --prefix      Prefix for outputs (e.g., AthCEN178)
  --consensus   Consensus FASTA used in Step 5 (comparemonomer)
  --seq         Genome FASTA (Step 2 arrayfind)
  --variant     Consensus variant FASTA (Step 1 getphmm)

Optional arguments:
  --bin         Bin size (default: 10000)
  --threads     Threads for consensusHORmonomer (default: 10)
  --plotv       Variant for sharedHOR plotting (default: V2)  [V1|V2|V3]
  --min-monomer-len  Minimum monomer length for monomerfind (default: 120)
  -h, --help    Show this help and exit

Notes:
  • Step 5 will ALWAYS include --consensus as requested.
  • Step 6 (classprediction) runs WITHOUT --plot by default.
  • This orchestrator assumes subcommands are already on PATH as 'HiReNET'.
============================================================
EOF
}

# ---------- defaults ----------
CHR_CSV=""
PREFIX=""
CONSENSUS=""
GENOME=""
VARIANT=""
BIN=10000
THREADS=10
PLOTV="V2"
MIN_MONOMER_LEN=120

# ---------- parse args ----------
[[ $# -eq 0 ]] && { usage; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chr)        CHR_CSV="${2:-}"; shift 2;;
    --prefix)     PREFIX="${2:-}"; shift 2;;
    --consensus)  CONSENSUS="${2:-}"; shift 2;;
    --seq)        GENOME="${2:-}"; shift 2;;
    --variant)    VARIANT="${2:-}"; shift 2;;
    --bin)        BIN="${2:-}"; shift 2;;
    --threads)    THREADS="${2:-}"; shift 2;;
    --plotv)      PLOTV="${2:-}"; shift 2;;
    --min-monomer-len) MIN_MONOMER_LEN="${2:-}"; shift 2;;
    -h|--help)    usage; exit 0;;
    *)            die "Unknown argument: $1 (see --help)";;
  esac
done

# ---------- validate ----------
[[ -n "$CHR_CSV"   ]] || die "Missing --chr"
[[ -n "$PREFIX"    ]] || die "Missing --prefix"
[[ -n "$CONSENSUS" ]] || die "Missing --consensus"
[[ -n "$GENOME"    ]] || die "Missing --seq"
[[ -n "$VARIANT"   ]] || die "Missing --variant"

# Normalize paths (requires readlink -f; matches your other scripts)
CONSENSUS="$(readlink -f "$CONSENSUS")"
GENOME="$(readlink -f "$GENOME")"
VARIANT="$(readlink -f "$VARIANT")"

# ---------- derived standard paths ----------
PHMM_DIR="phmm"
ARRAY_DIR="${PREFIX}_arrayout"
MONO_DIR="${PREFIX}_monomerout"
ARR_BIN_DIR="${PREFIX}_arrangemonomer_${BIN/000/}kb"
CMP1_DIR="${PREFIX}_comparemonomers"
CLSP_DIR="${PREFIX}_classpred_out"
REARR_DIR="${PREFIX}_rearrange_monomers_mergebin"
CMP2_DIR="${PREFIX}_compare_rearrangemonomers"
NET_DIR="${PREFIX}_network_HOR_mergebin"
ARR_HOR_DIR="${PREFIX}_network_mergebin_consensus"

# ---------- pipeline ----------
echo "=== Step 1 — getphmm ==="
#ml MAFFT/7.526-GCC-13.3.0-with-extensions
ml HMMER/3.4-gompi-2023b
#HiReNET getphmm -i "$VARIANT" -o "$PHMM_DIR" -p "$PREFIX"

echo "=== Step 2 — arrayfind ==="
ml BEDTools/2.31.1-GCC-13.3.0
ml BLAST+/2.16.0-gompi-2024a
#ml SeqKit/2.9.0
#HiReNET arrayfind -g "$GENOME" -c "$CONSENSUS" -o "$ARRAY_DIR" -p "$PREFIX"

echo "=== Step 3 — monomerfind ==="
ml HMMER/3.4-gompi-2023b
#ml BEDTools/2.31.1-GCC-13.3.0
ml bioawk/1.0-GCC-12.3.0
HiReNET monomerfind \
  --arrays-dir "$ARRAY_DIR" \
  --chrom-dir "$ARRAY_DIR/split_seq" \
  --outdir "$MONO_DIR" \
  --prefix "$PREFIX" \
  --hmm "${PREFIX}_$PHMM_DIR/${PREFIX}.hmm" \
  --chr "$CHR_CSV" \
  --min-monomer-len "$MIN_MONOMER_LEN"

echo "=== Step 4 — arrangemonomer ==="
HiReNET arrangemonomer \
  --arrays-dir "$ARRAY_DIR" \
  --genomic-bed-dir "$MONO_DIR" \
  --monomer-dir "$MONO_DIR/${PREFIX}_monomers" \
  --outdir "$ARR_BIN_DIR" \
  --prefix "$PREFIX" \
  --bin "$BIN" \
  --chr "$CHR_CSV"

echo "=== Step 5 — comparemonomer (self + consensus) ==="
ml BLAT/3.7-GCC-12.3.0
HiReNET comparemonomer \
  --bins-dir "$ARR_BIN_DIR/${PREFIX}_bin_monomers" \
  --outdir "$CMP1_DIR" \
  --consensus "$CONSENSUS"

echo "=== Step 6 — classprediction (no --plot) ==="
ml R/4.3.1-foss-2022a
HiReNET classprediction \
  --blatsub "$CMP1_DIR/blat_output_sub" \
  --outdir "$CLSP_DIR" \
  --prefix "$PREFIX" \
  --bin "$BIN"

echo "=== Step 7 — rearrangemonomers (HOR) ==="
HiReNET rearrangemonomers \
  --bins "$CLSP_DIR/${PREFIX}_fin_bins_combined.txt" \
  --class HOR \
  --prefix "$PREFIX" \
  --monomer-dir "$MONO_DIR/${PREFIX}_monomers" \
  --outdir "$REARR_DIR" \
  --chr "$CHR_CSV"

echo "=== Step 8 — comparemonomer on rearranged ==="
HiReNET comparemonomer \
  --bins-dir "$REARR_DIR/re_arrange_monomers" \
  --outdir "$CMP2_DIR"

echo "=== Step 9 — networkHOR (merged HOR bins) ==="
HiReNET networkHOR \
  --blatsub "$CMP2_DIR/blat_output_sub" \
  --bins "$CLSP_DIR/${PREFIX}_fin_bins_combined.txt" \
  --coor "$REARR_DIR/${PREFIX}_monomer_bed_inbin.txt" \
  --outdir "$NET_DIR"

echo "=== Step 10 — arrangeHORmonomer ==="
HiReNET arrangeHORmonomer \
  --groupdir "$NET_DIR" \
  --monomer-dir "$MONO_DIR/${PREFIX}_monomers" \
  --outdir "$ARR_HOR_DIR"

echo "=== Step 11 — per-chromosome consensus, compare, shared (plotv=${PLOTV}) ==="
IFS=',' read -ra CHRS <<< "$CHR_CSV"
for chr in "${CHRS[@]}"; do
  chr="$(echo "$chr" | xargs)"  # trim
  [[ -n "$chr" ]] || continue

  HiReNET consensusHORmonomer \
    --outdir "$ARR_HOR_DIR" \
    --threads "$THREADS" \
    --chroms "$chr"

  HiReNET compareConsensus \
    --chr "$chr" \
    --consensdir "$ARR_HOR_DIR/all_recluster_consensus_monomer" \
    --outdir "${PREFIX}_compare_consensusHOR_${chr}"

  HiReNET sharedHOR \
    --chr "$chr" \
    --datadir "${PREFIX}_compare_consensusHOR_${chr}/blat_sub" \
    --outdir "${PREFIX}_shared_out_${chr}" \
    --letter "${NET_DIR}/mergebin_string_outputs" \
    --plotv "$PLOTV"
done

echo "[✓] runALL completed."