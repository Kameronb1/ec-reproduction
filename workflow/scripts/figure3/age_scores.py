"""Figure 3H, step 2: the fetal-to-adult age score of every drawn row.

Each row is one profile: the pooled raw counts of its own cells, as log2(CPM + 1), with CPM over the 21,520 genes the
intestinal pseudobulk measures (every other input measures them too, so all profiles share one denominator). Each list
gene is centered on the mean of the four intestinal arm profiles, so gene baseline levels cancel and 0 lies midway
between fetal and adult intestinal EC. Score = mean centered value of the adult-high genes minus that of the fetal-high
genes, on the list with the 10% detection floor (7 adult-high, 9 fetal-high genes; age_genes.R).

Rows: four intestinal arms by primary study label (fetal EC of the Gut Cell Atlas pair, the mean of the Elmentaite 2020
and 2021 fetal profiles; Egozi 2021 fetal; Elmentaite 2021 adult; Hickey 2023), three intestinal test sources by
the four-marker classifier (Yu 2021, Wang 2020, Burclaff 2022: cells classified as EC; not used in the gene fit), SC-SIC and
SC-beta of each SC-islet study, Singh 2024 induced EC (day 6) and the McCauley 2023 HIO graft EC (Endocrine_1, both
organoids). Every profile needs >= 20 cells.

Inputs (name=path):
  gut_sums gut_groups     primary intestinal EEC pseudobulk by primary study label
  test_sums test_groups   intestinal test sources by four-marker classifier result (CALLED / ABSTAIN / NOT_CALLED)
  sc_sums sc_groups ind_sums ind_groups hio_sums hio_groups
  adult fetal             the gene lists with the detection floor
Outputs:
  drawn   the values panel H draws: block, row, kind of row, score, row position from the top
  genes   every drawn row gene by gene: own log2(CPM + 1), centered value and contribution to the score
"""
import sys

import numpy as np
import pandas as pd

arg = dict(a.split("=", 1) for a in sys.argv[1:])
MIN_CELLS = 20
LISTS = {"adult": open(arg["adult"]).read().split(), "fetal": open(arg["fetal"]).read().split()}
GENES = LISTS["adult"] + LISTS["fetal"]
E20, E21, HICKEY, EGOZI = ("SF_ELMENTAITE2020_DEVGUT", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC", "SF_EGOZI2021_DEVELOPMENTAL_EEC")
STUDY = {"GSE151117": "Augsornworawat 2020", "GSE167880": "Balboa 2022", "GSE185036": "Du 2022",
         "GSE199636": "Augsornworawat 2023 (nuclei)", "GSE202497": "Zhu 2023", "GSE226346": "Siehler 2024"}


def load(tag):
    s = pd.read_csv(arg[tag + "_sums"], sep="\t", index_col="gene")
    g = pd.read_csv(arg[tag + "_groups"], sep="\t", dtype={"unit": str})
    assert list(s.columns) == list(g.group) and not s.index.has_duplicates, tag
    return s, g.set_index("group")


PB = {t: load(t) for t in ("gut", "test", "sc", "ind", "hio")}
UNIVERSE = PB["gut"][0].index
assert all(UNIVERSE.isin(s.index).all() for s, _ in PB.values()), "an input misses a gene of the intestinal pseudobulk"
assert all(x in UNIVERSE for x in GENES), "a list gene is not measured"
if sorted(PB["sc"][1].family_id.unique()) != sorted(STUDY):
    sys.exit("the SC pseudobulk does not hold the six 10x studies")


def own(tag, mask):
    """log2(CPM + 1) of the list genes in the pooled counts of the groups in mask."""
    s, g = PB[tag]
    if g.loc[mask, "n_cells"].sum() < MIN_CELLS:
        sys.exit(f"a {tag} profile has fewer than {MIN_CELLS} cells")
    c = s.loc[UNIVERSE, g.index[mask]].sum(axis=1)
    return np.log2(c[GENES] / c.sum() * 1e6 + 1)


gg = PB["gut"][1]
arm = lambda fam, pred=lambda u: True: own("gut", (gg.family_id == fam) & (gg["class"] == "POSITIVE") & gg.unit.map(pred).astype(bool))
ARMS = {"Elmentaite 2020 + 2021 fetal": (arm(E20, lambda u: u[0].isdigit()) + arm(E21, lambda u: u.startswith("F"))) / 2,
        "Egozi 2021 fetal": arm(EGOZI, lambda u: u.endswith("w")),
        "Elmentaite 2021 adult": arm(E21, lambda u: u.startswith("A")),
        "Hickey 2023 adult (nuclei)": arm(HICKEY)}
CENTER = sum(ARMS.values()) / len(ARMS)

tg = PB["test"][1]
TEST = {"Yu 2021 fetal *": ("SF_YU2021_ENDODERM", "fetal gut EC"), "Wang 2020 adult *": ("SF_GSE125970_ADULT_INTESTINE", "adult gut EC"),
        "Burclaff 2022 adult *": ("SF_GSE185224_ADULT_REGIONAL_INTESTINE", "adult gut EC")}
sg, ig, hg = PB["sc"][1], PB["ind"][1], PB["hio"][1]
profiles = [("Gut tissue", r, "fetal gut EC" if "fetal" in r else "adult gut EC", v) for r, v in ARMS.items()]
profiles += [("Gut tissue", r, what, own("test", (tg.family_id == fam) & (tg["class"] == "CALLED"))) for r, (fam, what) in TEST.items()]
for pop, what in (("SC_SIC", "SC-SIC"), ("SC_BETA", "SC-beta")):
    profiles += [("Grown in vitro", STUDY[st], what, own("sc", (sg.family_id == st) & (sg["class"] == pop))) for st in STUDY]
profiles += [("Grown in vitro", "Singh 2024 induced EEC, day 6", "induced EEC", own("ind", (ig.unit == "D6") & (ig["class"] == "EC"))),
             (" ", "McCauley 2023 HIO graft EC", "HIO graft EC", own("hio", hg.unit.isin(["wt1", "wt2"]) & (hg["class"] == "EC")))]

genes, drawn = [], []
for block, row, what, level in profiles:
    centered = level - CENTER
    contrib = pd.concat([centered[LISTS["adult"]] / len(LISTS["adult"]), -centered[LISTS["fetal"]] / len(LISTS["fetal"])])
    genes += [{"row": row, "what": what, "gene": x, "list": "adult-high" if x in LISTS["adult"] else "fetal-high",
               "log2cpm": level[x], "centered": centered[x], "contribution": contrib[x]} for x in GENES]
    drawn.append({"block": block, "row": row, "what": what, "score": float(contrib.sum())})
D = pd.DataFrame(drawn)

# row order, top to bottom: intestinal rows adult first, each age by score (highest at the top); then the SC-islet studies by
# their SC-SIC score (a study's SC-SIC and SC-beta share a row); then induced EC and the HIO graft
gut = D[D.block == "Gut tissue"].sort_values(["what", "score"], ascending=[True, False])     # "adult ..." sorts before "fetal ..."
sic = D[D.what == "SC-SIC"].sort_values("score", ascending=False)
pos = {(r, w): i for i, (r, w) in enumerate(zip(gut.row, gut.what), 1)}
spos = {r: i for i, r in enumerate(sic.row, len(gut) + 1)}
rest = D[(D.block != "Gut tissue") & ~D.what.isin(["SC-SIC", "SC-beta"])]
pos.update({(r, w): i for i, (r, w) in enumerate(zip(rest.row, rest.what), len(gut) + len(sic) + 1)})
D["row_from_top"] = [spos[r] if w in ("SC-SIC", "SC-beta") else pos[(r, w)] for r, w in zip(D.row, D.what)]
D.to_csv(arg["drawn"], sep="\t", index=False, lineterminator="\n", float_format="%.15g")
pd.DataFrame(genes).to_csv(arg["genes"], sep="\t", index=False, lineterminator="\n", float_format="%.15g")
