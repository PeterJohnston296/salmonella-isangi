#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${PHENIX_INPUT_DIR}/phenix_reads.tsv}"
OUTPUT_DIR="${2:-${PHENIX_OUTPUT_DIR}}"
require_tsv_header "$MANIFEST" $'analysis_id\tr1\tr2'
require_file "$REFERENCE_FASTA"
require_file "$PHENIX_CONFIG"
require_command nextflow
require_command docker
mkdir -p "$OUTPUT_DIR"

COMMAND_LOG="${OUTPUT_DIR}/nextflow_command.log"
args=(nextflow run "${SCRIPT_DIR}/09_phenix.nf" -with-docker "$PHENIX_CONTAINER"
  --reads_manifest "$MANIFEST"
  --reference "$REFERENCE_FASTA"
  --phenix_config "$PHENIX_CONFIG"
  --output_dir "$OUTPUT_DIR"
  --container_image "$PHENIX_CONTAINER"
  --max_forks "$PHENIX_MAX_FORKS")
if [[ "${PHENIX_RESUME:-0}" == "1" ]]; then args+=(-resume); fi
run_logged "$COMMAND_LOG" "${args[@]}"
write_versions "${OUTPUT_DIR}/software_versions.tsv" nextflow docker
note "PHEnix run completed. Preserve the Nextflow report, trace and container digest before release."
