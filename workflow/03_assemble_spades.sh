#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${TRIM_DIR}/trimmed_reads.tsv}"
OUTPUT_DIR="${2:-${ASSEMBLY_DIR}}"
require_tsv_header "$MANIFEST" $'analysis_id\tr1\tr2'
require_command spades.py
mkdir -p "$OUTPUT_DIR"

OUTPUT_MANIFEST="${OUTPUT_DIR}/assemblies.tsv"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
assert_new_path "$OUTPUT_MANIFEST"
printf 'analysis_id\tassembly\n' > "$OUTPUT_MANIFEST"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r analysis_id r1 r2; do
  safe_analysis_id "$analysis_id"
  require_file "$r1"
  require_file "$r2"
  sample_dir="${OUTPUT_DIR}/${analysis_id}_assembly"
  contigs="${sample_dir}/contigs.fasta"
  assert_new_path "$sample_dir"
  run_logged "$COMMAND_LOG" spades.py -1 "$r1" -2 "$r2" \
    -o "$sample_dir" --careful -t "$SPADES_THREADS"
  require_file "$contigs"
  printf '%s\t%s\n' "$analysis_id" "$contigs" >> "$OUTPUT_MANIFEST"
done

write_versions "${OUTPUT_DIR}/software_versions.tsv" spades.py
note "Assembly manifest: $OUTPUT_MANIFEST"
