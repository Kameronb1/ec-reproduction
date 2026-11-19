"""Figure S9A (McCauley 2023 HIO graft EC column): check of the selected candidate regulators (Figure 5B) in McCauley
2023 human intestinal organoid grafts (GSE214852), on the published labels.

What it computes. For every regulon of every selected factor (program_values from regulon_selection.py): whether the
regulon is scored in GSE214852 (coverage >= 10 measured targets and >= 50% of its targets, network not trained on
GSE214852), and its effect = mean over the 3 graft preparations of (mean regulon AUC in EC cells minus mean in other
enteroendocrine cells), with EC = the published Endocrine_1 cluster (322 cells) and other enteroendocrine = Endocrine_2
(276 cells). The matched-null p value is not computed for the published labels (written NA). Not used in factor
selection.
Inputs (name=path): program_values, effects (coverage and network of each regulon in GSE214852), training, cells, auc
Output: out
"""
import csv, gzip, io, sys
from collections import defaultdict

a = dict(x.split("=", 1) for x in sys.argv[1:])
TARGET = "GSE214852__primary"


def read(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def tsv_text(rows):
    """The rows as a tab-separated table (csv.DictWriter), read back as text, so every value is the text of the published table."""
    h = io.StringIO()
    w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    return list(csv.DictReader(io.StringIO(h.getvalue()), delimiter="\t"))


# ---- step 1: the effect of every regulon from per-cell AUC with the published labels -------------------------------------
cells = {r["cell_id"]: r for r in read(a["cells"])}
label = {c: ("POSITIVE" if r["author_label"] == "Endocrine_1" else "CONTROL") for c, r in cells.items()}
assert {k: sum(v == k for v in label.values()) for k in ("POSITIVE", "CONTROL")} == {"POSITIVE": 322, "CONTROL": 276}

auc_by_net = defaultdict(list)   # network -> [(cell, regulon, auc text)] in file order
for r in read(a["auc"]):
    auc_by_net[r["network"]].append((r["cell_id"], r["regulon"], r["auc"]))


def effect(recs):
    acc = defaultdict(lambda: [0.0, 0])
    for c, p, v in recs:
        k = acc[(p, cells[c]["sample"], label[c])]
        k[0] += float(v); k[1] += 1
    means = {k: v[0] / v[1] for k, v in acc.items()}
    samples = sorted({k[1] for k in means})
    progs = list(dict.fromkeys(p for _, p, _ in recs))
    return {p: sum(means[(p, s, "POSITIVE")] - means[(p, s, "CONTROL")] for s in samples) / len(samples) for p in progs}


effects_primary = {}
for net, recs in auc_by_net.items():
    effects_primary.update(effect(recs))

# ---- step 2: which selected regulons are scored in GSE214852 (coverage, training) ---------------------------------------
training = {r["network"]: r["training_studies"].split(",") for r in read(a["training"])}
cov = {r["regulon"]: r for r in read(a["effects"]) if r["target"] == TARGET}
rows = []
for r in read(a["program_values"]):
    z = cov.get(r["program"])
    x = None if z is None else dict(cov=z["coverage_pass"] == "TRUE", external="GSE214852" not in training[z["network"]])
    ok = x is not None and x["cov"] and x["external"] and r["program"] in effects_primary
    rows.append(dict(factor=r["factor"], direction=r["direction"], d=r["d"], program=r["program"], network=r["network"],
                     program_pass=r["program_pass"], gse214852_scored=ok,
                     gse214852_effect=repr(effects_primary[r["program"]]) if ok else None,
                     gse214852_matched_p="NA" if ok else None,
                     gse214852_coverage_pass=bool(x and x["cov"]), gse214852_external=bool(x and x["external"])))
out_rows = tsv_text(rows)
with open(a["out"], "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, list(out_rows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(out_rows)
print("scored regulons:", sum(r["gse214852_scored"] == "True" for r in out_rows))
