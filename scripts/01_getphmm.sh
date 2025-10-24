#!/usr/bin/env bash
set -euo pipefail

# --- find the HiReNET folder, but don't change working directory ---
if [[ -n "${HIRENET_ROOT:-}" && -d "${HIRENET_ROOT}" ]]; then
  ROOT="${HIRENET_ROOT}"
else
  ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
export HIRENET_ROOT="${ROOT}"
# stay in the user’s working directory

# ==============================================================
# HiReNET: getphmm.sh  (Linux version)
# --------------------------------------------------------------
# Purpose:
#   Build a nucleotide profile Hidden Markov Model (HMM) from a
#   DNA FASTA sequence. This script aligns the input sequence(s)
#   using MAFFT, trims the alignment if necessary, and constructs
#   an HMM profile with HMMER.
#
# Dependencies:
#   - mafft
#   - hmmbuild (HMMER)
#   - seqkit (optional, for FASTA validation)
#
# Input:
#   - DNA FASTA file (single or multiple sequences)
#
# Output:
#   - <outdir>/<prefix>.aln.fasta   : MAFFT alignment
#   - <outdir>/<prefix>.hmm         : Profile HMM model
#
# Example:
#   HiReNET getphmm \
#       -i AthCEN178_consensus_variant.fasta \
#       -o AthCEN178_phmm \
#       -p AthCEN178 \
#       -t 8
#
# Notes:
#   - Automatically creates the output directory if missing.
#   - Detects available CPU threads if -t is not specified.
#   - Produces .aln.fasta suitable for downstream monomer detection.
# ==============================================================

IN_FASTA=""
OUTBASE=""
PREFIX=""
THREADS=""
DNA="--dna"   # set empty if protein

usage() {
  cat <<'EOF'
============================================================
                     HiReNET: getphmm
============================================================

Description:
  Build a profile HMM from a DNA FASTA sequence using MAFFT 
  for alignment and HMMER for model generation.

------------------------------------------------------------
Usage:
  HiReNET getphmm -i <input.fasta> -o <outdir> -p <prefix> [options]

Required arguments:
  -i   Input FASTA file (DNA)
  -o   Output directory
  -p   Prefix for output files

Optional arguments:
  -t   Number of threads for MAFFT (default: auto-detect)
  -h   Show this help message and exit

------------------------------------------------------------
Outputs:
  <outdir>/<prefix>.aln.fasta      Aligned sequences
  <outdir>/<prefix>.hmm            HMM profile

Example:
  HiReNET getphmm -i AthCEN178_consensus_variant.fasta \
                  -o AthCEN178_phmm \
                  -p AthCEN178

============================================================
EOF
}

# --- Parse arguments ---
while getopts "i:o:p:t:h" opt; do
  case "$opt" in
    i) IN_FASTA="$OPTARG" ;;
    o) OUTBASE="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# --- Validate ---
[[ -z "$IN_FASTA" || -z "$OUTBASE" || -z "$PREFIX" ]] && { echo "[x] Missing -i/-o/-p"; usage; }
[[ -f "$IN_FASTA" ]] || { echo "[x] Input FASTA not found: $IN_FASTA"; exit 1; }

# --- Tool checks (works on Mac/conda) ---
command -v mafft >/dev/null 2>&1 || { echo "[x] mafft not found in PATH"; exit 1; }
command -v hmmbuild >/dev/null 2>&1 || { echo "[x] hmmbuild not found in PATH"; exit 1; }

# --- Threads (auto for Mac) ---
if [[ -z "${THREADS}" ]]; then
  if command -v sysctl >/dev/null 2>&1; then
    THREADS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
  else
    THREADS=4
  fi
fi

# --- Output dir: PREFIX_OUTBASE ---
OUTDIR="${PREFIX}_${OUTBASE}"
mkdir -p "$OUTDIR"

ALN="${OUTDIR}/${PREFIX}.aln.fasta"
HMM_OUT="${OUTDIR}/${PREFIX}.hmm"

echo "[*] MAFFT threads: ${THREADS}"
echo "[*] Aligning ${IN_FASTA} → ${ALN}"
mafft --thread "${THREADS}" --auto "${IN_FASTA}" > "${ALN}"

echo "[*] Building HMM (${DNA}) → ${HMM_OUT}"
hmmbuild ${DNA} "${HMM_OUT}" "${ALN}"

echo "[✓] Done"
echo "    Output dir : ${OUTDIR}"
echo "    Alignment  : ${ALN}"
echo "    HMM        : ${HMM_OUT}"
