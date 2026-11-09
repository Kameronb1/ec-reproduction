"""Check every input listed in a figure's manifest against its recorded SHA-256.

Usage: verify_inputs.py <inputs/figure.tsv> <data_dir> <report.txt>
Exits non-zero, naming each missing or changed file, if anything differs.
"""
import csv
import hashlib
import sys
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(manifest, data_dir, report):
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    problems, lines = [], []
    for row in rows:
        path = Path(data_dir) / row["file"]
        if not path.is_file():
            problems.append(f"missing: {row['file']}")
            continue
        observed = sha256(path)
        if observed != row["sha256"]:
            problems.append(f"changed: {row['file']} expected {row['sha256']} observed {observed}")
        else:
            lines.append(f"ok\t{row['file']}\t{observed}")
    if problems:
        sys.exit("input check failed for " + manifest + "\n" + "\n".join(problems))
    Path(report).parent.mkdir(parents=True, exist_ok=True)
    Path(report).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
