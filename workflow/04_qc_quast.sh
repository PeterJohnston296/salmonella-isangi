#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${ASSEMBLIES_MANIFEST}}"
OUTPUT_DIR="${2:-${QC_DIR}}"
require_tsv_header "$MANIFEST" $'analysis_id\tassembly'
require_nonempty QUAST_LENGTH_METRIC "$QUAST_LENGTH_METRIC"
require_nonempty QUAST_CONTIG_METRIC "$QUAST_CONTIG_METRIC"
require_command quast
require_command python3
mkdir -p "$OUTPUT_DIR"

REPORT_MANIFEST="${OUTPUT_DIR}/quast_reports.tsv"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
assert_new_path "$REPORT_MANIFEST"
printf 'analysis_id\treport\n' > "$REPORT_MANIFEST"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r analysis_id assembly; do
  safe_analysis_id "$analysis_id"
  require_file "$assembly"
  sample_dir="${OUTPUT_DIR}/${analysis_id}"
  report="${sample_dir}/report.tsv"
  assert_new_path "$sample_dir"
  run_logged "$COMMAND_LOG" quast -o "$sample_dir" -t "$QUAST_THREADS" "$assembly"
  require_file "$report"
  printf '%s\t%s\n' "$analysis_id" "$report" >> "$REPORT_MANIFEST"
done

python3 "${SCRIPT_DIR}/lib/summarise_quast.py" \
  "$REPORT_MANIFEST" "${OUTPUT_DIR}/assembly_qc_audit.tsv" \
  --length-metric "$QUAST_LENGTH_METRIC" \
  --contig-metric "$QUAST_CONTIG_METRIC" \
  --max-length "$MAX_ASSEMBLY_LENGTH" \
  --max-contigs "$MAX_CONTIGS"

write_versions "${OUTPUT_DIR}/software_versions.tsv" quast
note "QC audit written without changing the frozen membership lists: ${OUTPUT_DIR}/assembly_qc_audit.tsv"
