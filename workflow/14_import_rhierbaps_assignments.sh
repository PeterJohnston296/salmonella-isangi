#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

ASSIGNMENTS="${1:-}"
MEMBERSHIP="${2:-${CC25_IDS}}"
OUTPUT_DIR="${3:-${PHYLOGENY_DIR}/rhierbaps}"
require_file "$ASSIGNMENTS"
require_file "$MEMBERSHIP"
require_command python3
assert_new_path "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
cp "$ASSIGNMENTS" "${OUTPUT_DIR}/rhierbaps_assignments.tsv"
python3 "${SCRIPT_DIR}/lib/validate_assignments.py" \
  "${OUTPUT_DIR}/rhierbaps_assignments.tsv" "$MEMBERSHIP" "${OUTPUT_DIR}/validation.json"
note "Imported frozen lineage assignments. RHierBAPS was not rerun by this repository."
