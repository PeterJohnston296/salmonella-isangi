#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

ALIGNMENT="${1:-}"
OUTPUT="${2:-}"
require_file "$ALIGNMENT"
require_nonempty OUTPUT "$OUTPUT"
require_command snp-dists
assert_new_path "$OUTPUT"
mkdir -p "$(dirname "$OUTPUT")"

record_command "${OUTPUT}.command.log" snp-dists -j "$SNP_DISTS_THREADS" "$ALIGNMENT"
snp-dists -j "$SNP_DISTS_THREADS" "$ALIGNMENT" > "$OUTPUT"
{
  printf 'field\tvalue\n'
  printf 'alignment\t%s\n' "$ALIGNMENT"
  printf 'alignment_sha256\t%s\n' "$(sha256_file "$ALIGNMENT")"
  printf 'matrix\t%s\n' "$OUTPUT"
  printf 'matrix_sha256\t%s\n' "$(sha256_file "$OUTPUT")"
  printf 'note\tSNP distances were calculated directly from the named input alignment; recombination status is not inferred from its filename.\n'
} > "${OUTPUT}.provenance.tsv"
note "SNP matrix: $OUTPUT"
