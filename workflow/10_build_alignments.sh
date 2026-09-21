#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

SUBSET="${1:-}"
MEMBERSHIP="${2:-}"
PHENIX_ROOT="${3:-${PHENIX_OUTPUT_DIR}}"
OUTPUT_DIR="${4:-${ALIGNMENT_DIR}/${SUBSET}}"
[[ "$SUBSET" =~ ^[A-Za-z0-9_.-]+$ ]] || die "Usage: $0 SUBSET MEMBERSHIP_FILE [PHENIX_ROOT] [OUTPUT_DIR]"
require_file "$MEMBERSHIP"
require_file "$REFERENCE_HEADERS_FILE"
require_command python3
assert_new_path "$OUTPUT_DIR"
mkdir -p "${OUTPUT_DIR}/sample_only"

SAMPLE_MANIFEST="${OUTPUT_DIR}/sample_fastas.tsv"
printf 'analysis_id\tfasta\n' > "$SAMPLE_MANIFEST"
while IFS= read -r analysis_id; do
  [[ -n "$analysis_id" ]] || continue
  safe_analysis_id "$analysis_id"
  source_fasta="${PHENIX_ROOT}/${analysis_id}/phenix_bbduk/${analysis_id}_all.fasta"
  require_file "$source_fasta"
  sample_fasta="${OUTPUT_DIR}/sample_only/${analysis_id}.fasta"
  python3 "${SCRIPT_DIR}/lib/fasta_tools.py" extract \
    "$source_fasta" "$analysis_id" "$REFERENCE_HEADERS_FILE" "$sample_fasta" \
    >> "${OUTPUT_DIR}/extraction_records.jsonl"
  printf '%s\t%s\n' "$analysis_id" "$sample_fasta" >> "$SAMPLE_MANIFEST"
done < "$MEMBERSHIP"

ALIGNMENT="${OUTPUT_DIR}/${SUBSET}_reference_coordinate_alignment.fasta"
python3 "${SCRIPT_DIR}/lib/fasta_tools.py" combine "$SAMPLE_MANIFEST" "$MEMBERSHIP" "$ALIGNMENT" \
  > "${OUTPUT_DIR}/combination_record.json"
python3 "${SCRIPT_DIR}/lib/fasta_tools.py" check "$ALIGNMENT" "$MEMBERSHIP" \
  "${OUTPUT_DIR}/alignment_check.json"
note "Alignment created with exact membership and row order: $ALIGNMENT"
