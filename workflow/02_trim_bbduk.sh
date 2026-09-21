#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${READS_MANIFEST}}"
OUTPUT_DIR="${2:-${TRIM_DIR}}"
require_tsv_header "$MANIFEST" $'analysis_id\tr1\tr2'
require_file "$ADAPTERS_FASTA"
require_command bbduk.sh
mkdir -p "$OUTPUT_DIR"

OUTPUT_MANIFEST="${OUTPUT_DIR}/trimmed_reads.tsv"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
assert_new_path "$OUTPUT_MANIFEST"
printf 'analysis_id\tr1\tr2\n' > "$OUTPUT_MANIFEST"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r analysis_id r1 r2; do
  safe_analysis_id "$analysis_id"
  require_file "$r1"
  require_file "$r2"
  out1="${OUTPUT_DIR}/${analysis_id}_1_trimmed.fastq.gz"
  out2="${OUTPUT_DIR}/${analysis_id}_2_trimmed.fastq.gz"
  assert_new_path "$out1"
  assert_new_path "$out2"
  run_logged "$COMMAND_LOG" bbduk.sh \
    "threads=${BBDUK_THREADS}" "ref=${ADAPTERS_FASTA}" \
    "in=${r1}" "in2=${r2}" "out=${out1}" "out2=${out2}" \
    ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=20 minlength=50
  [[ -s "$out1" && -s "$out2" ]] || die "BBDuk did not produce both outputs for $analysis_id"
  printf '%s\t%s\t%s\n' "$analysis_id" "$out1" "$out2" >> "$OUTPUT_MANIFEST"
done

write_versions "${OUTPUT_DIR}/software_versions.tsv" bbduk.sh
note "Trimmed-read manifest: $OUTPUT_MANIFEST"
