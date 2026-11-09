"""Download the inputs that come from their original public source rather than the deposit.

Usage (from the repository root): python workflow/scripts/common/fetch_url_inputs.py <data_dir>
Reads every inputs/*.tsv row with origin "url", downloads it into <data_dir>/<file> unless a copy
with the right checksum is already there, and refuses a download whose SHA-256 differs.
"""
import csv
import glob
import hashlib
import sys
import urllib.request
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(data_dir):
    wanted = {}
    for manifest in sorted(glob.glob("inputs/*.tsv")):
        with open(manifest, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["origin"] == "url":
                    wanted[row["file"]] = row
    for name, row in wanted.items():
        target = Path(data_dir) / name
        if target.is_file() and sha256(target) == row["sha256"]:
            print(f"present: {name}")
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_name(target.name + ".part")
        print(f"downloading {row['url']}")
        urllib.request.urlretrieve(row["url"], partial)
        if sha256(partial) != row["sha256"]:
            partial.unlink()
            sys.exit(f"{name}: downloaded file does not match the recorded SHA-256; nothing kept")
        partial.replace(target)
        print(f"ok: {name}")


if __name__ == "__main__":
    main(sys.argv[1])
