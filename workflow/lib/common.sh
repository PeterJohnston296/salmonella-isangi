#!/usr/bin/env bash


die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

note() {
  printf '[isangi] %s\n' "$*" >&2
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "Install the executable for this stage: $1"
}

require_file() {
  [[ -f "$1" ]] || die "Stage input: $1"
}

require_directory() {
  [[ -d "$1" ]] || die "Stage directory: $1"
}

require_nonempty() {
  [[ -n "${2:-}" ]] || die "Required setting is empty: $1"
}

safe_analysis_id() {
  [[ "$1" =~ ^[A-Za-z0-9_.#-]+$ ]] || die "Unsafe or unsupported analysis ID: $1"
}

assert_new_path() {
  local target="$1"
  if [[ -e "$target" && "${ALLOW_EXISTING_OUTPUTS:-0}" != "1" ]]; then
    die "Output already exists: $target. Set ALLOW_EXISTING_OUTPUTS=1 only after checking it."
  fi
}

require_tsv_header() {
  local file="$1"
  local expected="$2"
  require_file "$file"
  local actual
  IFS= read -r actual < "$file" || die "Could not read header from $file"
  actual="${actual%$'\r'}"
  [[ "$actual" == "$expected" ]] || die "Unexpected header in $file. Expected: $expected. Found: $actual"
}

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  else
    die "Neither sha256sum nor shasum is available"
  fi
}

record_command() {
  local log_file="$1"
  shift
  mkdir -p "$(dirname "$log_file")"
  {
    printf '%q ' "$@"
    printf '\n'
  } >> "$log_file"
}

run_logged() {
  local log_file="$1"
  shift
  record_command "$log_file" "$@"
  "$@"
}

write_versions() {
  local output="$1"
  shift
  mkdir -p "$(dirname "$output")"
  : > "$output"
  local tool
  for tool in "$@"; do
    if command -v "$tool" >/dev/null 2>&1; then
      {
        printf '%s\t' "$tool"
        "$tool" --version 2>&1 | head -n 1 || true
      } >> "$output"
    else
      printf '%s\tnot_installed\n' "$tool" >> "$output"
    fi
  done
}
