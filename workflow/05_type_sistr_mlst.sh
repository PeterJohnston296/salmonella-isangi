#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${ASSEMBLIES_MANIFEST}}"
OUTPUT_DIR="${2:-${TYPING_DIR}}"
require_tsv_header "$MANIFEST" $'analysis_id\tassembly'
require_command sistr
require_command mlst
require_command python3
mkdir -p "$OUTPUT_DIR"

SISTR_MANIFEST="${OUTPUT_DIR}/sistr_reports.tsv"
MLST_COMBINED="${OUTPUT_DIR}/mlst_results.tsv"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
assert_new_path "$SISTR_MANIFEST"
assert_new_path "$MLST_COMBINED"
printf 'analysis_id\treport\n' > "$SISTR_MANIFEST"
printf 'analysis_id\tmlst_output\n' > "$MLST_COMBINED"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r analysis_id assembly; do
  safe_analysis_id "$analysis_id"
  require_file "$assembly"
  sample_dir="${OUTPUT_DIR}/${analysis_id}"
  mkdir -p "$sample_dir"
  sistr_output="${sample_dir}/sistr.csv"
  mlst_output="${sample_dir}/mlst.tsv"
  assert_new_path "$sistr_output"
  assert_new_path "$mlst_output"
  run_logged "$COMMAND_LOG" sistr -i "$assembly" "$analysis_id" -f csv -o "$sistr_output"
  record_command "$COMMAND_LOG" mlst --quiet "$assembly"
  mlst --quiet "$assembly" > "$mlst_output"
  printf '%s\t%s\n' "$analysis_id" "$sistr_output" >> "$SISTR_MANIFEST"
  mlst_line="$(tr '\t' ' ' < "$mlst_output" | tr -s ' ')"
  printf '%s\t%s\n' "$analysis_id" "$mlst_line" >> "$MLST_COMBINED"
done

python3 "${SCRIPT_DIR}/lib/merge_reports.py" \
  "$SISTR_MANIFEST" "${OUTPUT_DIR}/sistr_combined.csv" "${OUTPUT_DIR}/sistr_completion.tsv" \
  --delimiter comma
write_versions "${OUTPUT_DIR}/software_versions.tsv" sistr mlst
note "Typing results: $OUTPUT_DIR"
