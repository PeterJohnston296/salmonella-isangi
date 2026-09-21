#!/usr/bin/env python3
"""Validate a frozen one-row-per-isolate assignment table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

def read_ids(path: Path) -> list[str]:
    values = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    if not values or len(values) != len(set(values)):
        raise ValueError(f"Membership file is empty or contains duplicates: {path}")
    return values

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assignments", type=Path)
    parser.add_argument("membership", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()
    ids = read_ids(args.membership)
    with args.assignments.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["analysis_id", "lineage"]:
            parser.error("Assignments must contain exactly: analysis_id, lineage")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        parser.error("Ragged assignment table")
    row_ids = [row["analysis_id"] for row in rows]
    if len(row_ids) != len(set(row_ids)):
        parser.error("Duplicate analysis IDs in assignments")
    if row_ids != ids:
        parser.error("Assignment membership/order does not exactly match the supplied ID list")
    if any(not row["lineage"].strip() for row in rows):
        parser.error("Blank lineage assignment")
    result = {
        "assignments": str(args.assignments),
        "membership": str(args.membership),
        "n_isolates": len(rows),
        "n_lineages": len({row["lineage"] for row in rows}),
        "ordered_membership_verified": True,
        "assignments_sha256": hashlib.sha256(args.assignments.read_bytes()).hexdigest(),
        "membership_sha256": hashlib.sha256(args.membership.read_bytes()).hexdigest(),
        "method_status": "imported_frozen_assignments; clustering was not rerun",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
