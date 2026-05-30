#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build"
OUT_DIR="${SCRIPT_DIR}/output"
BIN="${BUILD_DIR}/easymsa-prep"
WORK_DIR="${OUT_DIR}/work"

echo "[1/8] Configure"
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}"

echo "[2/8] Build"
cmake --build "${BUILD_DIR}" -j

echo "[3/8] Prepare example archives"
rm -rf "${OUT_DIR}"
mkdir -p "${WORK_DIR}/archive/a" "${WORK_DIR}/archive/b"
cp "${SCRIPT_DIR}/mixed_sequences.fasta" "${WORK_DIR}/mixed_sequences.fasta"
cp "${SCRIPT_DIR}/mixed_sequences.fasta" "${WORK_DIR}/archive/a/part1.fasta"
cp "${SCRIPT_DIR}/rna_sequences.fasta" "${WORK_DIR}/archive/b/part2.fasta"
xz -zkf "${WORK_DIR}/mixed_sequences.fasta"
tar -czf "${WORK_DIR}/mixed_sequences.tar.gz" -C "${WORK_DIR}/archive" .
tar -cJf "${WORK_DIR}/mixed_sequences.tar.xz" -C "${WORK_DIR}/archive" .

cat > "${WORK_DIR}/rna_params.json" <<'JSON'
{
  "mode": "audit",
  "strictness": "normal",
  "sequence_type": "rna",
  "rules": {
    "u_to_t": true
  }
}
JSON

echo "[4/8] Run gzip input"
"${BIN}" "${SCRIPT_DIR}/mixed_sequences.fasta.gz" \
  -o "${OUT_DIR}/gz" \
  --mode audit \
  --strictness normal \
  --threads 4 \
  --force

echo "[5/8] Run xz input"
"${BIN}" "${WORK_DIR}/mixed_sequences.fasta.xz" \
  -o "${OUT_DIR}/xz" \
  --mode audit \
  --strictness normal \
  --threads 4 \
  --force

echo "[6/8] Run archive inputs"
"${BIN}" "${SCRIPT_DIR}/mixed_sequences.zip" \
  -o "${OUT_DIR}/zip" \
  --mode audit \
  --strictness normal \
  --threads 4 \
  --force

"${BIN}" "${WORK_DIR}/mixed_sequences.tar.gz" \
  -o "${OUT_DIR}/tar_gz" \
  --mode audit \
  --strictness normal \
  --threads 4 \
  --force

"${BIN}" "${WORK_DIR}/mixed_sequences.tar.xz" \
  -o "${OUT_DIR}/tar_xz" \
  --mode audit \
  --strictness normal \
  --threads 4 \
  --force

echo "[7/8] Run RNA params and deterministic thread comparison"
"${BIN}" "${SCRIPT_DIR}/rna_sequences.fasta" \
  -o "${OUT_DIR}/rna" \
  --params "${WORK_DIR}/rna_params.json" \
  --threads 4 \
  --force

"${BIN}" "${SCRIPT_DIR}/mixed_sequences.fasta" \
  -o "${OUT_DIR}/t1" \
  --threads 1 \
  --force

"${BIN}" "${SCRIPT_DIR}/mixed_sequences.fasta" \
  -o "${OUT_DIR}/t8" \
  --threads 8 \
  --force

cmp "${OUT_DIR}/t1/clean.fasta" "${OUT_DIR}/t8/clean.fasta"

echo "[8/8] Show representative outputs"
echo
echo "=== Gzip summary ==="
cat "${OUT_DIR}/gz/summary.txt"

echo
echo "=== Gzip result JSON ==="
cat "${OUT_DIR}/gz/result.json"

echo
echo "=== Clean FASTA ==="
cat "${OUT_DIR}/gz/clean.fasta"

echo
echo "=== Output tree ==="
find "${OUT_DIR}" -maxdepth 3 -type f | sort

echo
echo "Examples completed successfully."
