#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"
ALIGNMENT="${1:-}"
NAME="${2:-isangi}"
OUTPUT_DIR="${3:-${PHYLOGENY_DIR}/classic_raxml}"
require_file "$ALIGNMENT"
require_command raxmlHPC-PTHREADS
ALIGNMENT="$(realpath "$ALIGNMENT")"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"
run_logged "${OUTPUT_DIR}/commands.log" raxmlHPC-PTHREADS -f a -m GTRGAMMA \
  -p 12345 -x 12345 -# 1000 -T "$RAXML_THREADS" -s "$ALIGNMENT" -n "$NAME" -w "$OUTPUT_DIR"
write_versions "${OUTPUT_DIR}/software_versions.tsv" raxmlHPC-PTHREADS
