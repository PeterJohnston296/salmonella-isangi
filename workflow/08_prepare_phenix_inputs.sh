#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

TRIMMED_MANIFEST="${1:-${TRIM_DIR}/trimmed_reads.tsv}"
MEMBERSHIP="${2:-${PHENIX_COHORT_IDS}}"
OUTPUT_DIR="${3:-${PHENIX_INPUT_DIR}}"
require_tsv_header "$TRIMMED_MANIFEST" $'analysis_id\tr1\tr2'
require_file "$MEMBERSHIP"
require_file "$REFERENCE_FASTA"
require_file "$PHENIX_CONFIG"
assert_new_path "$OUTPUT_DIR"
mkdir -p "${OUTPUT_DIR}/reads"

declare -A R1_BY_ID R2_BY_ID
while IFS=$'\t' read -r analysis_id r1 r2; do
  [[ "$analysis_id" == "analysis_id" ]] && continue
  safe_analysis_id "$analysis_id"
  [[ -z "${R1_BY_ID[$analysis_id]+x}" ]] || die "Duplicate read ID: $analysis_id"
  R1_BY_ID[$analysis_id]="$r1"
  R2_BY_ID[$analysis_id]="$r2"
done < "$TRIMMED_MANIFEST"

OUTPUT_MANIFEST="${OUTPUT_DIR}/phenix_reads.tsv"
printf 'analysis_id\tr1\tr2\n' > "$OUTPUT_MANIFEST"
while IFS= read -r analysis_id; do
  [[ -n "$analysis_id" ]] || continue
  safe_analysis_id "$analysis_id"
  r1="${R1_BY_ID[$analysis_id]:-}"
  r2="${R2_BY_ID[$analysis_id]:-}"
  [[ -n "$r1" && -n "$r2" ]] || die "PHEnix cohort ID is absent from trimmed reads: $analysis_id"
  require_file "$r1"
  require_file "$r2"
  link1="${OUTPUT_DIR}/reads/${analysis_id}_1_trimmed.fastq.gz"
  link2="${OUTPUT_DIR}/reads/${analysis_id}_2_trimmed.fastq.gz"
  assert_new_path "$link1"
  assert_new_path "$link2"
  ln -s "$(realpath "$r1")" "$link1"
  ln -s "$(realpath "$r2")" "$link2"
  printf '%s\t%s\t%s\n' "$analysis_id" "$link1" "$link2" >> "$OUTPUT_MANIFEST"
done < "$MEMBERSHIP"

{
  printf 'resource\tpath\tsha256\n'
  printf 'reference\t%s\t%s\n' "$REFERENCE_FASTA" "$(sha256_file "$REFERENCE_FASTA")"
  printf 'phenix_config\t%s\t%s\n' "$PHENIX_CONFIG" "$(sha256_file "$PHENIX_CONFIG")"
  printf 'membership\t%s\t%s\n' "$MEMBERSHIP" "$(sha256_file "$MEMBERSHIP")"
} > "${OUTPUT_DIR}/input_checksums.tsv"
note "PHEnix input manifest: $OUTPUT_MANIFEST"
