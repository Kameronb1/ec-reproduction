"""Rebuild PANELS.tsv, the one-row-per-panel index, from panels/<figure>.tsv.

Usage (from the repository root): python workflow/scripts/common/build_panels_index.py
Each row follows one panel from its inputs to its image: inputs manifest, rules file, analysis
scripts, tables (and how many are checked against expected/<figure>.tsv), plot script, drawn image
and published image. It stops if a plot script or published image named in panels/<figure>.tsv
is missing.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
UNITS = ["figure1", "figure2", "figure2_supp", "figure3", "figure3_supp", "figure4", "figure5"]


def rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


out = []
for unit in UNITS:
    expected = {r["output"] for r in rows(ROOT / "expected" / f"{unit}.tsv")}
    for r in rows(ROOT / "panels" / f"{unit}.tsv"):
        outputs = [o.strip() for o in r["outputs"].split(",") if o.strip()]
        checked = [o for o in outputs if o in expected]
        plot = f"workflow/scripts/{unit}/{r['plot']}" if r["plot"] else ""
        published = f"expected/{unit}/panels/{r['image']}" if r["image"] else ""
        for path in (plot, published):
            if path and not (ROOT / path).is_file():
                raise SystemExit(f"{r['panel']}: {path} is missing")
        out.append({
            "panel": r["panel"],
            "shows": r["shows"],
            "inputs": f"inputs/{unit}.tsv",
            "rules": f"workflow/rules/{unit}.smk",
            "scripts": r["script"],
            "tables": ",".join(f"results/{unit}/{o}" for o in outputs),
            "checked_against_published": f"{len(checked)} of {len(outputs)}" if outputs else "no table (a key)",
            "plot": plot,
            "image": f"results/{unit}/panels/{r['image']}" if r["image"] else "",
            "published_image": published,
        })

with open(ROOT / "PANELS.tsv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(out)
print(f"PANELS.tsv: {len(out)} panels")
