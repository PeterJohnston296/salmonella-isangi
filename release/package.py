#!/usr/bin/env python3
from pathlib import Path
import csv
import hashlib
import zipfile

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {".git", ".Rproj.user", ".quarto", "__pycache__", "work", "library", "staging"}
PUBLIC_DIRS = {"data", "scripts", "workflow", "provenance", "docs", "tests", "environment", ".github", "release"}
PUBLIC_FILES = {"README.md", "NOTICE.md", "LICENSE", "CITATION.cff", ".zenodo.json", ".gitignore", ".gitattributes", "Makefile", "_quarto.yml", "salmonella-isangi.Rproj"}


def include(path):
    relative = path.relative_to(ROOT)
    if EXCLUDED.intersection(relative.parts):
        return False
    if path.name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".cram", ".ttf", ".otf", ".woff", ".woff2", ".partial")):
        return False
    return (
        relative.parts[0] in PUBLIC_DIRS
        or relative.as_posix() in PUBLIC_FILES
        or relative.as_posix().startswith(("results/reference/", "results/tables/"))
    )


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


checksum = ROOT / "provenance/SHA256SUMS.txt"
file_index = ROOT / "docs/file_index.csv"

# Refresh the human-readable inventory before calculating release checksums.
indexed = sorted(
    path for path in ROOT.rglob("*")
    if path.is_file() and include(path) and path not in {checksum, file_index}
)
with file_index.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(("path", "bytes"))
    writer.writerows((path.relative_to(ROOT).as_posix(), path.stat().st_size) for path in indexed)

files = sorted(path for path in ROOT.rglob("*") if path.is_file() and include(path))
checksum.write_text(
    "".join(
        f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}\n"
        for path in files if path != checksum
    ),
    encoding="utf-8",
)
if checksum not in files:
    files.append(checksum)

archive = ROOT.parent / "salmonella-isangi_v1.0.0.zip"
with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
    for path in sorted(files):
        bundle.write(path, Path(ROOT.name) / path.relative_to(ROOT))

sidecar = archive.with_suffix(archive.suffix + ".sha256")
sidecar.write_text(f"{sha256(archive)}  {archive.name}\n", encoding="utf-8")
print(archive)
print(sidecar)
