#!/usr/bin/env python3
"""Combine one-assembly QUAST reports and audit manuscript thresholds."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

def digest(path: Path) -> str:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return h

def read_report(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or any(len(row) != 2 for row in rows):
        raise ValueError(f"Expected a single-assembly two-column QUAST report: {path}")
    metrics = {row[0]: row[1] for row in rows[1:]}
    if len(metrics) != len(rows) - 1:
        raise ValueError(f"Duplicate QUAST metric names in {path}")
    return metrics

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV: analysis_id, report")
    parser.add_argument("output", type=Path)
    parser.add_argument("--length-metric", required=True)
    parser.add_argument("--contig-metric", required=True)
    parser.add_argument("--max-length", type=int, required=True)
    parser.add_argument("--max-contigs", type=int, required=True)
    args = parser.parse_args()
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["analysis_id", "report"]:
            parser.error("Manifest must contain exactly: analysis_id, report")
        rows = list(reader)
    output_rows = []
    for row in rows:
        path = Path(row["report"])
        metrics = read_report(path)
        try:
            length = int(metrics[args.length_metric])
            contigs = int(metrics[args.contig_metric])
        except KeyError as error:
            parser.error(f"Metric {error.args[0]!r} is absent from {path}")
        output_rows.append({
            "analysis_id": row["analysis_id"],
            "length_metric": args.length_metric,
            "length_bp": length,
            "contig_metric": args.contig_metric,
            "contigs": contigs,
            "exceeds_manuscript_limits": str(length > args.max_length or contigs > args.max_contigs).lower(),
            "report_sha256": digest(path),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        fields = list(output_rows[0]) if output_rows else ["analysis_id"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
