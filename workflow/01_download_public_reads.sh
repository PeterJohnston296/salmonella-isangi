#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${1:-${REPO_ROOT}/data/manifests/public_archive_runs_274.notebook_derived.tsv}"
OUTPUT_DIR="${2:-${RAW_READS_DIR}}"
require_tsv_header "$MANIFEST" $'sample_accession\trun_accession\tevidence_status'
require_command fastq-dump
require_command gzip

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
RESOLVED="${OUTPUT_DIR}/resolved_runs.tsv"
COMMAND_LOG="${OUTPUT_DIR}/commands.log"
assert_new_path "$RESOLVED"
printf 'sample_accession\trun_accession\tr1\tr2\n' > "$RESOLVED"

resolve_runs() {
  local sample_accession="$1"
  require_command esearch
  require_command efetch
  esearch -db sra -query "$sample_accession" \
    | efetch -format runinfo \
    | cut -d, -f1 \
    | grep -E '^(SRR|ERR|DRR)[0-9]+$' \
    | sort -u
}

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample_accession run_accession evidence_status; do
  [[ -n "$sample_accession" || -n "$run_accession" ]] || die "Manifest row has neither sample nor run accession"
  runs="$run_accession"
  if [[ -z "$runs" ]]; then
    runs="$(resolve_runs "$sample_accession")"
    [[ -n "$runs" ]] || die "No run accession resolved for $sample_accession"
  fi
  while IFS= read -r run; do
    [[ "$run" =~ ^(SRR|ERR|DRR)[0-9]+$ ]] || die "Unexpected run accession: $run"
    r1="${OUTPUT_DIR}/${run}_1.fastq.gz"
    r2="${OUTPUT_DIR}/${run}_2.fastq.gz"
    if [[ -f "$r1" && -f "$r2" ]]; then
      note "Retaining existing read pair for $run"
    else
      assert_new_path "$r1"
      assert_new_path "$r2"
      (
        cd "$OUTPUT_DIR"
        run_logged "$COMMAND_LOG" fastq-dump --split-files "$run"
        [[ -s "${run}_1.fastq" && -s "${run}_2.fastq" ]] || die "Expected paired FASTQ files were not produced for $run"
        run_logged "$COMMAND_LOG" gzip "${run}_1.fastq" "${run}_2.fastq"
      )
    fi
    printf '%s\t%s\t%s\t%s\n' "$sample_accession" "$run" "$r1" "$r2" >> "$RESOLVED"
  done <<< "$runs"
done

write_versions "${OUTPUT_DIR}/software_versions.tsv" fastq-dump esearch efetch
note "Resolved-read manifest: $RESOLVED"
