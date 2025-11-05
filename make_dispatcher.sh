#!/usr/bin/env bash
set -euo pipefail

# Define main folders relative to this repo
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${ROOT}/bin"
SCRIPTS="${ROOT}/scripts"
RDIR="${ROOT}/R"

mkdir -p "${BIN}"

cat > "${BIN}/HiReNET" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

# --- find project root ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SCRIPTS="${ROOT}/scripts"
RDIR="${ROOT}/R"
export HIRENET_ROOT="${ROOT}"

usage() {
  echo "HiReNET — run HOR and monomer analysis pipeline"
  echo
  echo "Usage:"
  echo "  HiReNET <command> [args...]"
  echo
  echo "Available bash commands:"
  for f in $(ls -1 "${SCRIPTS}"/[0-9][0-9]_*.sh 2>/dev/null | sort); do
    b=$(basename "$f")
    printf "  %-22s (%s)\n" "$(sed -E 's/^[0-9]+_//' <<< "${b%.sh}")" "$b"
  done
  echo
  echo "Available R commands:"
  for f in $(ls -1 "${RDIR}"/S*_*.R 2>/dev/null | sort); do
    b=$(basename "$f")
    printf "  %-22s (%s)\n" "$(sed -E 's/^S[0-9]+_//' <<< "${b%.R}")" "$b"
  done
  echo
  echo "Example:"
  echo "  HiReNET getphmm -i input.fasta -o out -p tag"
}

cmd="${1:-}"
[[ -z "$cmd" || "$cmd" == "-h" || "$cmd" == "--help" ]] && { usage; exit 0; }
shift || true

# Try Bash scripts first
file_bash=$(ls "${SCRIPTS}"/[0-9][0-9]_"${cmd}".sh 2>/dev/null | head -n1 || true)
if [[ -n "$file_bash" ]]; then
  exec bash "$file_bash" "$@"
fi

# Try R scripts
file_r=$(ls "${RDIR}"/S*_"${cmd}".R 2>/dev/null | head -n1 || true)
if [[ -n "$file_r" ]]; then
  exec Rscript "$file_r" "$@"
fi

echo "Unknown command: $cmd"
echo
usage
exit 1
EOF

chmod +x "${BIN}/HiReNET"
echo "Dispatcher created at ${BIN}/HiReNET"
echo "Add to PATH with (bash):"
echo "    echo 'export PATH=\"${BIN}:\$PATH\"' >> ~/.bashrc && source ~/.bashrc"
echo "Or for zsh:"
echo "    echo 'export PATH=\"${BIN}:\$PATH\"' >> ~/.zshrc && source ~/.zshrc"