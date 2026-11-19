"""Rebuild PANELS.tsv, the one-row-per-panel index, from panels/<figure>.tsv.

Usage (from the repository root): python workflow/scripts/common/build_panels_index.py
Each row follows one panel from its inputs to its image: inputs manifest, rules file, analysis
scripts, tables (and how many are checked against the expected tables), plot script, drawn image
and published image. A table or script written as ../<figure>/<name> in panels/<figure>.tsv belongs
to that other figure: its table is made by that figure's rules and checked against that figure's
expected tables. The script stops if a plot script, analysis script or published image named in
panels/<figure>.tsv is missing.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
UNITS = ["figure1", "figure2", "figure3", "figure4", "figure5", "figure_s1", "figure_s2_s3", "figure_s4",
         "figure_s5", "figure_s6", "figure_s7", "figure_s8", "figure_s9"]


def rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def owner(unit, name):
    """(figure, name) of a table or script; ../<figure>/<name> names another figure's."""
    if name.startswith("../"):
        other, rest = name[3:].split("/", 1)
        return other, rest
    return unit, name


def names(cell):
    return [o.strip() for o in cell.split(",") if o.strip()]


EXPECTED = {unit: {r["output"] for r in rows(ROOT / "expected" / f"{unit}.tsv")} for unit in UNITS}

out = []
for unit in UNITS:
    for r in rows(ROOT / "panels" / f"{unit}.tsv"):
        outputs = [owner(unit, o) for o in names(r["outputs"])]
        scripts = [owner(unit, s) for s in names(r["script"])]
        checked = [o for o in outputs if o[1] in EXPECTED[o[0]]]
        plot = f"workflow/scripts/{unit}/{r['plot']}" if r["plot"] else ""
        published = f"expected/{unit}/panels/{r['image']}" if r["image"] else ""
        script_paths = [f"workflow/scripts/{u}/{s}" for u, s in scripts]
        for path in [plot, published] + script_paths:
            if path and not (ROOT / path).is_file():
                raise SystemExit(f"{r['panel']}: {path} is missing")
        rule_units = [unit] + [u for u, _ in outputs + scripts if u != unit]
        out.append({
            "panel": r["panel"],
            "shows": r["shows"],
            "inputs": ",".join(f"inputs/{u}.tsv" for u in dict.fromkeys(rule_units)),
            "rules": ",".join(f"workflow/rules/{u}.smk" for u in dict.fromkeys(rule_units)),
            "scripts": ",".join(script_paths),
            "tables": ",".join(f"results/{u}/{o}" for u, o in outputs),
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
