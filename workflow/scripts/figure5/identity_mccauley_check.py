"""Figure S18 B (McCauley 2023 column; also the factor directions read by S18 A and S19): outside check of the selected
identity factors in McCauley 2023 human intestinal organoid grafts (GSE214852), on the primary study labels.

What it computes. For every regulon of every selected factor (program_values from identity_selection.py): whether the
regulon is scored in GSE214852 (coverage >= 10 measured targets and >= 50% of its targets, network not trained on
GSE214852), and its effect = mean over the 3 graft preparations of (mean regulon AUC in EC cells minus mean in other
enteroendocrine cells), with EC = the primary study's Endocrine_1 cluster (322 cells) and other enteroendocrine = Endocrine_2
(276 cells). The effect with EC defined as the cells classified as EC by the Figure 2 label-free four-marker EC
classifier (319 vs 279 cells) is also written, as gse214852_effect_marker_rule. The matched-null p value is not
recomputed for the primary labels (written NA). Not used in factor selection.
Inputs (name=path): program_values, effects, training, cells, auc
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


def fnum(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def tsv_text(rows):
    """The rows as a tab-separated table (csv.DictWriter), read back as text, so every value is the text of the published table."""
    h = io.StringIO()
    w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    return list(csv.DictReader(io.StringIO(h.getvalue()), delimiter="\t"))


# ---- step 1: which regulons are scored in GSE214852 (coverage, training) and their effect under the four-marker classifier labels ----
training = {r["network"]: r["training_studies"].split(",") for r in read(a["training"])}
eff = {}
for r in read(a["effects"]):
    if r["target"] == TARGET:
        eff[r["regulon"]] = r
pv = read(a["program_values"])
rows = []
for r in pv:
    z = eff.get(r["program"])
    x = None if z is None else dict(effect=fnum(z["effect"]), p=fnum(z["matched_p"]), cov=z["coverage_pass"] == "TRUE",
                                    external="GSE214852" not in training[z["network"]])
    ok = x is not None and x["cov"] and x["external"] and x["effect"] is not None
    rows.append(dict(factor=r["factor"], direction=r["direction"], d=r["d"], program=r["program"], network=r["network"],
                     program_pass=r["program_pass"], gse214852_scored=ok,
                     gse214852_effect=x["effect"] if ok else None, gse214852_matched_p=x["p"] if ok else None,
                     gse214852_coverage_pass=bool(x and x["cov"]), gse214852_external=bool(x and x["external"])))
scored_rows = tsv_text(rows)

# ---- step 2: recompute the effect from per-cell AUC with the primary study labels ----------------------------------------
# The two labelings differ in 3 cells, all EC under the primary labels and other enteroendocrine under the classifier.
cells = {r["cell_id"]: r for r in read(a["cells"])}
by_classifier = {c: ("POSITIVE" if r["population"] == "EC" else "CONTROL") for c, r in cells.items()}
by_primary = {c: ("POSITIVE" if r["author_label"] == "Endocrine_1" else "CONTROL") for c, r in cells.items()}
moved = [c for c in cells if by_classifier[c] != by_primary[c]]
assert len(moved) == 3 and all(by_classifier[c] == "CONTROL" and by_primary[c] == "POSITIVE" for c in moved)
assert {k: sum(v == k for v in by_primary.values()) for k in ("POSITIVE", "CONTROL")} == {"POSITIVE": 322, "CONTROL": 276}

auc_by_net = defaultdict(list)   # network -> [(cell, regulon, auc text)] in file order
for r in read(a["auc"]):
    auc_by_net[r["network"]].append((r["cell_id"], r["regulon"], r["auc"]))
effects_classifier, effects_primary = {}, {}


def effect(recs, lab):
    acc = defaultdict(lambda: [0.0, 0])
    for c, p, v in recs:
        k = acc[(p, cells[c]["sample"], lab[c])]
        k[0] += float(v); k[1] += 1
    means = {k: v[0] / v[1] for k, v in acc.items()}
    samples = sorted({k[1] for k in means})
    progs = list(dict.fromkeys(p for _, p, _ in recs))
    return {p: sum(means[(p, s, "POSITIVE")] - means[(p, s, "CONTROL")] for s in samples) / len(samples) for p in progs}


for net, recs in auc_by_net.items():
    e_cls, e_pri = effect(recs, by_classifier), effect(recs, by_primary)
    for p in e_cls:
        if p in eff:
            effects_classifier[p] = e_cls[p]; effects_primary[p] = e_pri[p]

out_rows = []
for r in scored_rows:
    r = dict(r)
    if r["gse214852_scored"] == "True":
        p = r["program"]
        assert abs(effects_classifier[p] - float(r["gse214852_effect"])) <= 1e-12   # reproduces the stored classifier-label effect
        r["gse214852_effect_marker_rule"] = r["gse214852_effect"]
        r["gse214852_effect"] = repr(effects_primary[p])
        r["gse214852_matched_p"] = "NA"
    else:
        r["gse214852_effect_marker_rule"] = r["gse214852_effect"]
    out_rows.append(r)
cols = list(scored_rows[0].keys()) + ["gse214852_effect_marker_rule"]
with open(a["out"], "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(out_rows)
print("scored regulons:", sum(r["gse214852_scored"] == "True" for r in out_rows))
