#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# getphmm.sh
# Build profile HMM from input FASTA.
#
# Usage:
#   ./getphmm.sh -i input.fasta -o output_dir -p prefix [-t threads]
#
# Example:
#   ./getphmm.sh -i AthCEN178_consensus_variant.fasta -o phmm -p AthCEN178
#   -> Outputs in: AthCEN178_phmm/
# =========================================================

IN_FASTA=""
OUTBASE=""
PREFIX=""
THREADS=16

usage() {
  cat <<EOF
Usage: $0 -i input.fasta -o output_dir -p prefix [-t threads] [-h]

Required arguments:
  -i    Input FASTA file
  -o    Output directory base (prefix will be prepended)
  -p    Prefix name

Optional arguments:
  -t    Number of threads for MAFFT (default: 16)
  -h    Show this help message and exit

Outputs:
  PREFIX_OUTBASE/PREFIX.aln.fasta
  PREFIX_OUTBASE/PREFIX.hmm
EOF
  exit 0
}

# --- Parse arguments ---
while getopts "i:o:p:t:h" opt; do
  case "$opt" in
    i) IN_FASTA="$OPTARG" ;;
    o) OUTBASE="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

if [[ -z "${IN_FASTA}" || -z "${OUTBASE}" || -z "${PREFIX}" ]]; then
  echo "Error: missing arguments."
  usage
fi

# --- Load modules (adjust/remove to match your env) ---
module load MAFFT/7.526-GCC-13.3.0-with-extensions
module load HMMER/3.4-gompi-2023a

# --- Setup output directory with prefix ---
OUTDIR="${OUTBASE}"
mkdir -p "$OUTDIR"

ALN="${OUTDIR}/${PREFIX}.aln.fasta"
HMM_OUT="${OUTDIR}/${PREFIX}.hmm"

# --- Align ---
echo "Running MAFFT on ${IN_FASTA} ..."
mafft --thread "${THREADS}" "${IN_FASTA}" > "${ALN}"

# --- Build HMM ---
echo "Building HMM model ..."
hmmbuild "${HMM_OUT}" "${ALN}"

echo "Done!"
echo "  Output dir : ${OUTDIR}"
echo "  Alignment  : ${ALN}"
echo "  HMM        : ${HMM_OUT}"