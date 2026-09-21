#!/usr/bin/env python3
"""Strict FASTA extraction, combination and membership checks.

This release utility replaces the historical substring-based awk filter with
exact header matching. It does not infer aliases or silently drop records.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

class InputError(ValueError):
    pass

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def fasta_records(path: Path):
    name = None
    sequence = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if not sequence:
                        raise InputError(f"Empty sequence {name!r} in {path}")
                    yield name, "".join(sequence)
                name = line[1:].split()[0]
                if not name:
                    raise InputError(f"Blank FASTA header in {path}")
                sequence = []
            else:
                if name is None:
                    raise InputError(f"Sequence before first header in {path}")
                if not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn?\-]+", line):
                    raise InputError(f"Non-nucleotide characters in {path}, record {name}")
                sequence.append(line)
    if name is None:
        raise InputError(f"No FASTA records in {path}")
    if not sequence:
        raise InputError(f"Empty final sequence {name!r} in {path}")
    yield name, "".join(sequence)

def write_record(handle, name: str, sequence: str) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")

def read_ids(path: Path) -> list[str]:
    ids = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines()
           if line.strip() and not line.startswith("#")]
    if not ids:
        raise InputError(f"No IDs in {path}")
    if len(ids) != len(set(ids)):
        raise InputError(f"Duplicate IDs in {path}")
    return ids

def read_reference_headers(path: Path) -> set[str]:
    return set(read_ids(path))

def extract(source: Path, analysis_id: str, reference_headers: Path, output: Path) -> dict:
    records = list(fasta_records(source))
    names = [name for name, _ in records]
    if len(names) != len(set(names)):
        raise InputError(f"Duplicate FASTA headers in {source}")
    refs = read_reference_headers(reference_headers)
    if analysis_id in refs:
        raise InputError(f"Analysis ID {analysis_id} is also listed as a reference header")
    unexpected = [name for name in names if name != analysis_id and name not in refs]
    kept = [(name, seq) for name, seq in records if name == analysis_id]
    if unexpected:
        raise InputError(f"Unexpected FASTA records in {source}: {unexpected}")
    if len(kept) != 1:
        raise InputError(f"Expected exactly one record named {analysis_id!r} in {source}; found {names}")
    if len({len(seq) for _, seq in records}) != 1:
        raise InputError(f"Unequal record lengths in {source}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".partial")
    with temporary.open("w", encoding="utf-8") as handle:
        write_record(handle, *kept[0])
    temporary.replace(output)
    return {"source": str(source), "output": str(output), "analysis_id": analysis_id,
            "alignment_length": len(kept[0][1]), "sha256": sha256(output)}

def combine(manifest: Path, ids_file: Path, output: Path) -> dict:
    import csv

    ids = read_ids(ids_file)
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["analysis_id", "fasta"]:
            raise InputError(f"Expected analysis_id and fasta columns in {manifest}")
        rows = list(reader)
    lookup = {row["analysis_id"]: Path(row["fasta"]) for row in rows}
    if len(lookup) != len(rows):
        raise InputError(f"Duplicate IDs in {manifest}")
    if set(lookup) != set(ids):
        raise InputError(f"Manifest and membership list differ: missing={sorted(set(ids)-set(lookup))[:10]}, "
                         f"extra={sorted(set(lookup)-set(ids))[:10]}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".partial")
    alignment_length = None
    with temporary.open("w", encoding="utf-8") as handle:
        for analysis_id in ids:
            records = list(fasta_records(lookup[analysis_id]))
            if len(records) != 1 or records[0][0] != analysis_id:
                raise InputError(f"Expected one exact record {analysis_id!r} in {lookup[analysis_id]}")
            if alignment_length is None:
                alignment_length = len(records[0][1])
            if len(records[0][1]) != alignment_length:
                raise InputError(f"Alignment length differs for {analysis_id}")
            write_record(handle, *records[0])
    temporary.replace(output)
    return {"output": str(output), "n_sequences": len(ids),
            "alignment_length": alignment_length, "sha256": sha256(output),
            "ordered_ids_sha256": sha256(ids_file)}

def check(alignment: Path, ids_file: Path) -> dict:
    ids = read_ids(ids_file)
    records = list(fasta_records(alignment))
    names = [name for name, _ in records]
    if names != ids:
        raise InputError("FASTA row order or membership does not exactly match the supplied ID list")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise InputError("Alignment sequences have unequal lengths")
    return {"alignment": str(alignment), "n_sequences": len(records),
            "alignment_length": next(iter(lengths)), "ordered_membership_verified": True,
            "sha256": sha256(alignment), "ordered_ids_sha256": sha256(ids_file)}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("extract")
    p.add_argument("source", type=Path)
    p.add_argument("analysis_id")
    p.add_argument("reference_headers", type=Path)
    p.add_argument("output", type=Path)
    p = sub.add_parser("combine")
    p.add_argument("manifest", type=Path)
    p.add_argument("ids", type=Path)
    p.add_argument("output", type=Path)
    p = sub.add_parser("check")
    p.add_argument("alignment", type=Path)
    p.add_argument("ids", type=Path)
    p.add_argument("output_json", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "extract":
            result = extract(args.source, args.analysis_id, args.reference_headers, args.output)
        elif args.command == "combine":
            result = combine(args.manifest, args.ids, args.output)
        else:
            result = check(args.alignment, args.ids)
        if args.command == "check":
            args.output_json.parent.mkdir(parents=True, exist_ok=True)
            args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        else:
            print(json.dumps(result, sort_keys=True))
    except (InputError, OSError) as error:
        parser.exit(2, f"INPUT ERROR: {error}\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
