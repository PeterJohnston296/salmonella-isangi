#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

SUBSET="${1:-}"
ALIGNMENT="${2:-}"
OUTPUT_DIR="${3:-${PHYLOGENY_DIR}/${SUBSET}/raxmlng}"
[[ "$SUBSET" =~ ^[A-Za-z0-9_.-]+$ ]] || die "Usage: $0 SUBSET GUBBINS_FILTERED_ALIGNMENT [OUTPUT_DIR]"
require_file "$ALIGNMENT"
require_command raxml-ng
[[ "$TREE_MODEL" == "GTR+G" ]] || die "Manuscript wrapper requires TREE_MODEL=GTR+G"
[[ "$TREE_BOOTSTRAPS" == "1000" ]] || die "Manuscript wrapper requires TREE_BOOTSTRAPS=1000"
assert_new_path "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

PREFIX="${OUTPUT_DIR}/${SUBSET}"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
run_logged "$COMMAND_LOG" raxml-ng --all \
  --msa "$ALIGNMENT" \
  --model "$TREE_MODEL" \
  --prefix "$PREFIX" \
  --seed "$TREE_SEED" \
  --threads "$RAXML_THREADS" \
  --bs-trees "$TREE_BOOTSTRAPS"
require_file "${PREFIX}.raxml.bestTree"
require_file "${PREFIX}.raxml.support"
{
  printf 'field\tvalue\n'
  printf 'input_alignment\t%s\n' "$ALIGNMENT"
  printf 'input_sha256\t%s\n' "$(sha256_file "$ALIGNMENT")"
  printf 'best_tree\t%s\n' "${PREFIX}.raxml.bestTree"
  printf 'best_tree_sha256\t%s\n' "$(sha256_file "${PREFIX}.raxml.bestTree")"
  printf 'support_tree\t%s\n' "${PREFIX}.raxml.support"
  printf 'support_tree_sha256\t%s\n' "$(sha256_file "${PREFIX}.raxml.support")"
  printf 'model\t%s\n' "$TREE_MODEL"
  printf 'bootstrap_replicates\t%s\n' "$TREE_BOOTSTRAPS"
  printf 'seed\t%s\n' "$TREE_SEED"
  printf 'historical_equivalence\tnot established; late notebook commands used classic raxmlHPC-PTHREADS\n'
} > "${OUTPUT_DIR}/raxmlng_provenance.tsv"
write_versions "${OUTPUT_DIR}/software_versions.tsv" raxml-ng
note "RAxML-NG best and support trees: $OUTPUT_DIR"
