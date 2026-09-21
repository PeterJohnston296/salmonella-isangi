#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${ASSEMBLIES_MANIFEST}}"
MEMBERSHIP="${2:-${QC_APPROVED_IDS}}"
OUTPUT_DIR="${3:-${AMR_DIR}}"
require_tsv_header "$MANIFEST" $'analysis_id\tassembly'
require_file "$MEMBERSHIP"
require_directory "$AMRFINDER_DB"
require_command amrfinder
require_command python3
mkdir -p "$OUTPUT_DIR"

declare -A ASSEMBLY_BY_ID
while IFS=$'\t' read -r analysis_id assembly; do
  [[ "$analysis_id" == "analysis_id" ]] && continue
  safe_analysis_id "$analysis_id"
  [[ -z "${ASSEMBLY_BY_ID[$analysis_id]+x}" ]] || die "Duplicate assembly ID: $analysis_id"
  ASSEMBLY_BY_ID[$analysis_id]="$assembly"
done < "$MANIFEST"

REPORT_MANIFEST="${OUTPUT_DIR}/amrfinder_reports.tsv"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
assert_new_path "$REPORT_MANIFEST"
printf 'analysis_id\treport\n' > "$REPORT_MANIFEST"

while IFS= read -r analysis_id; do
  [[ -n "$analysis_id" ]] || continue
  safe_analysis_id "$analysis_id"
  assembly="${ASSEMBLY_BY_ID[$analysis_id]:-}"
  [[ -n "$assembly" ]] || die "No assembly for membership ID: $analysis_id"
  require_file "$assembly"
  sample_dir="${OUTPUT_DIR}/${analysis_id}"
  mkdir -p "$sample_dir"
  output="${sample_dir}/amrfinder.tsv"
  assert_new_path "$output"
  args=(amrfinder -n "$assembly" --organism Salmonella -d "$AMRFINDER_DB" -o "$output")
  if [[ "$AMRFINDER_USE_PLUS" == "1" ]]; then
    args+=(--plus)
  fi
  run_logged "$COMMAND_LOG" "${args[@]}"
  require_file "$output"
  printf '%s\t%s\n' "$analysis_id" "$output" >> "$REPORT_MANIFEST"
done < "$MEMBERSHIP"

python3 "${SCRIPT_DIR}/lib/merge_reports.py" \
  "$REPORT_MANIFEST" "${OUTPUT_DIR}/amrfinder_combined.tsv" "${OUTPUT_DIR}/amrfinder_completion.tsv" \
  --delimiter tab
write_versions "${OUTPUT_DIR}/software_versions.tsv" amrfinder
note "AMRFinder results retain explicit zero-hit completion records: $OUTPUT_DIR"
