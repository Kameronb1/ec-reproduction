"""Download the inputs that come from their original public source rather than the deposit.

Usage (from the repository root): python workflow/scripts/common/fetch_url_inputs.py <data_dir>
Reads every inputs/*.tsv row with origin "url", downloads it into <data_dir>/<file> unless a copy
with the right checksum is already there, and refuses a download whose SHA-256 differs.
"""
import csv
import glob
import hashlib
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, path, attempts=5):
    """Fetch url into path. Some servers refuse the default Python client name or limit the request
    rate, so the request names this workflow and a refused request is repeated after a pause."""
    request = urllib.request.Request(url, headers={"User-Agent": "ec-reproduction-workflow/1.0", "Accept": "*/*"})
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response, open(path, "wb") as handle:
                shutil.copyfileobj(response, handle)
            return
        except urllib.error.HTTPError as error:
            if error.code not in (429, 503) or attempt == attempts:
                sys.exit(f"{url}: the server answered {error.code}. Download the file by hand from this "
                         f"address, put it at the path given in inputs/, and run this script again.")
            wait = max(int(error.headers.get("Retry-After") or 0), 10 * attempt)
            print(f"  server answered {error.code}; trying again in {wait} s")
            time.sleep(wait)


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
        download(row["url"], partial)
        if sha256(partial) != row["sha256"]:
            partial.unlink()
            sys.exit(f"{name}: downloaded file does not match the recorded SHA-256; nothing kept")
        partial.replace(target)
        print(f"ok: {name}")


if __name__ == "__main__":
    main(sys.argv[1])
