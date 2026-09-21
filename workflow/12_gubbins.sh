#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

SUBSET="${1:-}"
ALIGNMENT="${2:-}"
OUTPUT_DIR="${3:-${PHYLOGENY_DIR}/${SUBSET}/gubbins}"
[[ "$SUBSET" =~ ^[A-Za-z0-9_.-]+$ ]] || die "Usage: $0 SUBSET ALIGNMENT [OUTPUT_DIR]"
require_file "$ALIGNMENT"
require_command run_gubbins.py
assert_new_path "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

PREFIX="${OUTPUT_DIR}/${SUBSET}_gubbins"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
run_logged "$COMMAND_LOG" run_gubbins.py --prefix "$PREFIX" --threads "$GUBBINS_THREADS" "$ALIGNMENT"
FILTERED="${PREFIX}.filtered_polymorphic_sites.fasta"
require_file "$FILTERED"
{
  printf 'field\tvalue\n'
  printf 'input_alignment\t%s\n' "$ALIGNMENT"
  printf 'input_sha256\t%s\n' "$(sha256_file "$ALIGNMENT")"
  printf 'filtered_alignment\t%s\n' "$FILTERED"
  printf 'filtered_sha256\t%s\n' "$(sha256_file "$FILTERED")"
} > "${OUTPUT_DIR}/gubbins_provenance.tsv"
write_versions "${OUTPUT_DIR}/software_versions.tsv" run_gubbins.py
note "Gubbins filtered alignment: $FILTERED"
