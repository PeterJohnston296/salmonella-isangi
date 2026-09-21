#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${ASSEMBLIES_MANIFEST}}"
MEMBERSHIP="${2:-${QC_APPROVED_IDS}}"
OUTPUT_DIR="${3:-${PLASMID_DIR}/mob_typer}"
require_tsv_header "$MANIFEST" $'analysis_id\tassembly'
require_file "$MEMBERSHIP"
require_command mob_typer
mkdir -p "$OUTPUT_DIR"

declare -A ASSEMBLY_BY_ID
while IFS=$'\t' read -r analysis_id assembly; do
  [[ "$analysis_id" == "analysis_id" ]] && continue
  safe_analysis_id "$analysis_id"
  ASSEMBLY_BY_ID[$analysis_id]="$assembly"
done < "$MANIFEST"

REPORT_MANIFEST="${OUTPUT_DIR}/mob_typer_reports.tsv"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
assert_new_path "$REPORT_MANIFEST"
printf 'analysis_id\treport\n' > "$REPORT_MANIFEST"
while IFS= read -r analysis_id; do
  [[ -n "$analysis_id" ]] || continue
  assembly="${ASSEMBLY_BY_ID[$analysis_id]:-}"
  [[ -n "$assembly" ]] || die "No assembly for $analysis_id"
  require_file "$assembly"
  sample_dir="${OUTPUT_DIR}/${analysis_id}"
  mkdir -p "$sample_dir"
  output="${sample_dir}/mob_typer_results.txt"
  assert_new_path "$output"
  run_logged "$COMMAND_LOG" mob_typer --infile "$assembly" --out_file "$output"
  require_file "$output"
  printf '%s\t%s\n' "$analysis_id" "$output" >> "$REPORT_MANIFEST"
done < "$MEMBERSHIP"
write_versions "${OUTPUT_DIR}/software_versions.tsv" mob_typer
note "MOB-typer per-isolate reports: $OUTPUT_DIR"
