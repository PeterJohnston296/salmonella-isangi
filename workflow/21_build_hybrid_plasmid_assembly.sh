#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

ANALYSIS_ID="${1:-}"
LONG_READS="${2:-}"
SHORT_R1="${3:-}"
SHORT_R2="${4:-}"
OUTPUT_DIR="${5:-${PLASMID_DIR}/${ANALYSIS_ID}}"
safe_analysis_id "$ANALYSIS_ID"
require_file "$LONG_READS"
require_file "$SHORT_R1"
require_file "$SHORT_R2"
require_command flye
require_command medaka_consensus
require_command pypolca
assert_new_path "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

COMMAND_LOG="${OUTPUT_DIR}/commands.log"
run_logged "$COMMAND_LOG" flye --nano-raw "$LONG_READS" --out-dir "${OUTPUT_DIR}/01_flye" \
  --genome-size 5m --threads "$REFERENCE_THREADS"
require_file "${OUTPUT_DIR}/01_flye/assembly.fasta"

medaka_args=(-i "$LONG_READS" -d "${OUTPUT_DIR}/01_flye/assembly.fasta" -o "${OUTPUT_DIR}/02_medaka_1" -t "$REFERENCE_THREADS")
if [[ -n "${MEDAKA_MODEL:-}" ]]; then medaka_args+=(-m "$MEDAKA_MODEL"); fi
run_logged "$COMMAND_LOG" medaka_consensus "${medaka_args[@]}"
medaka_args=(-i "$LONG_READS" -d "${OUTPUT_DIR}/02_medaka_1/consensus.fasta" -o "${OUTPUT_DIR}/03_medaka_2" -t "$REFERENCE_THREADS")
if [[ -n "${MEDAKA_MODEL:-}" ]]; then medaka_args+=(-m "$MEDAKA_MODEL"); fi
run_logged "$COMMAND_LOG" medaka_consensus "${medaka_args[@]}"

run_logged "$COMMAND_LOG" pypolca run -a "${OUTPUT_DIR}/03_medaka_2/consensus.fasta" \
  -1 "$SHORT_R1" -2 "$SHORT_R2" -t "$PYPOLCA_THREADS" -o "${OUTPUT_DIR}/04_pypolca_1" --careful
run_logged "$COMMAND_LOG" pypolca run -a "${OUTPUT_DIR}/04_pypolca_1/pypolca_corrected.fasta" \
  -1 "$SHORT_R1" -2 "$SHORT_R2" -t "$PYPOLCA_THREADS" -o "${OUTPUT_DIR}/05_pypolca_2" --careful
require_file "${OUTPUT_DIR}/05_pypolca_2/pypolca_corrected.fasta"
cp "${OUTPUT_DIR}/05_pypolca_2/pypolca_corrected.fasta" "${OUTPUT_DIR}/${ANALYSIS_ID}_hybrid_polished.fasta"
write_versions "${OUTPUT_DIR}/software_versions.tsv" flye medaka_consensus pypolca
note "Hybrid assembly created. Plasmid-contig selection, Circlator and pling remain separate unrecovered steps."
