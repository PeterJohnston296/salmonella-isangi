#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

LONG_READS="${1:-}"
SHORT_R1="${2:-}"
SHORT_R2="${3:-}"
OUTPUT_DIR="${4:-${REFERENCE_DIR}}"
require_file "$LONG_READS"
require_file "$SHORT_R1"
require_file "$SHORT_R2"
require_command flye
require_command medaka_consensus
require_command pypolca
assert_new_path "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

COMMAND_LOG="${OUTPUT_DIR}/commands.log"
FLYE_DIR="${OUTPUT_DIR}/01_flye"
MEDAKA1_DIR="${OUTPUT_DIR}/02_medaka_round_1"
MEDAKA2_DIR="${OUTPUT_DIR}/03_medaka_round_2"
POLCA1_DIR="${OUTPUT_DIR}/04_pypolca_round_1"
POLCA2_DIR="${OUTPUT_DIR}/05_pypolca_round_2"
FINAL="${OUTPUT_DIR}/CHF10J1_final_polished_reference.fasta"

run_logged "$COMMAND_LOG" flye --nano-raw "$LONG_READS" --out-dir "$FLYE_DIR" \
  --genome-size 5m --threads "$REFERENCE_THREADS"
require_file "${FLYE_DIR}/assembly.fasta"

medaka_args=(-i "$LONG_READS" -d "${FLYE_DIR}/assembly.fasta" -o "$MEDAKA1_DIR" -t "$REFERENCE_THREADS")
if [[ -n "${MEDAKA_MODEL:-}" ]]; then medaka_args+=(-m "$MEDAKA_MODEL"); fi
run_logged "$COMMAND_LOG" medaka_consensus "${medaka_args[@]}"
require_file "${MEDAKA1_DIR}/consensus.fasta"

medaka_args=(-i "$LONG_READS" -d "${MEDAKA1_DIR}/consensus.fasta" -o "$MEDAKA2_DIR" -t "$REFERENCE_THREADS")
if [[ -n "${MEDAKA_MODEL:-}" ]]; then medaka_args+=(-m "$MEDAKA_MODEL"); fi
run_logged "$COMMAND_LOG" medaka_consensus "${medaka_args[@]}"
require_file "${MEDAKA2_DIR}/consensus.fasta"

run_logged "$COMMAND_LOG" pypolca run -a "${MEDAKA2_DIR}/consensus.fasta" \
  -1 "$SHORT_R1" -2 "$SHORT_R2" -t "$PYPOLCA_THREADS" -o "$POLCA1_DIR" --careful
require_file "${POLCA1_DIR}/pypolca_corrected.fasta"

run_logged "$COMMAND_LOG" pypolca run -a "${POLCA1_DIR}/pypolca_corrected.fasta" \
  -1 "$SHORT_R1" -2 "$SHORT_R2" -t "$PYPOLCA_THREADS" -o "$POLCA2_DIR" --careful
require_file "${POLCA2_DIR}/pypolca_corrected.fasta"
cp "${POLCA2_DIR}/pypolca_corrected.fasta" "$FINAL"

{
  printf 'file\tsha256\n'
  printf '%s\t%s\n' "$FINAL" "$(sha256_file "$FINAL")"
} > "${OUTPUT_DIR}/reference_manifest.tsv"
write_versions "${OUTPUT_DIR}/software_versions.tsv" flye medaka_consensus pypolca
note "Reference built as $FINAL"

