#!/usr/bin/env python3
"""Combine per-isolate CSV/TSV reports while preserving zero-hit isolates."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV: analysis_id, report")
    parser.add_argument("output", type=Path)
    parser.add_argument("status", type=Path)
    parser.add_argument("--delimiter", choices=["tab", "comma"], default="tab")
    args = parser.parse_args()
    delimiter = "\t" if args.delimiter == "tab" else ","
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        manifest = csv.DictReader(handle, delimiter="\t")
        if manifest.fieldnames != ["analysis_id", "report"]:
            parser.error("Manifest must contain exactly: analysis_id, report")
        items = list(manifest)
    header = None
    combined = []
    status_rows = []
    seen = set()
    for item in items:
        analysis_id = item["analysis_id"]
        if analysis_id in seen:
            parser.error(f"Duplicate analysis ID: {analysis_id}")
        seen.add(analysis_id)
        path = Path(item["report"])
        if not path.is_file():
            parser.error(f"Missing report is not a zero-hit report: {path}")
        if path.stat().st_size == 0:
            fields = None
            rows = []
        else:
            with path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle, delimiter=delimiter)
                fields = reader.fieldnames
                rows = list(reader)
            if not fields:
                parser.error(f"Nonempty report has no parseable header: {path}")
            if header is None:
                header = fields
            elif fields != header:
                parser.error(f"Header differs in {path}")
            combined.extend({"analysis_id": analysis_id, **row} for row in rows)
        status_rows.append({
            "analysis_id": analysis_id,
            "status": "complete_with_hits" if rows else "complete_zero_hits",
            "n_report_rows": len(rows),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["analysis_id"] + (header or []),
                                delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(combined)
    with args.status.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["analysis_id", "status", "n_report_rows", "sha256"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(status_rows)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
