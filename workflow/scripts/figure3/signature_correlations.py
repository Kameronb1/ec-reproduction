"""Figure 3C: how closely each population's EC-versus-rest signature agrees with intestinal EC.

For every population (intestinal EC in each intestinal study, SC-SIC / SC-beta / SC-alpha in each SC-islet study,
induced EEC day 6, HIO graft EEC) a signature is the log2 fold change of pooled raw counts
(log2(CPM + 1), whole-library CPM) of that population against the rest of its sample. Each signature is
compared with the intestinal EC signature of each intestinal study by Spearman correlation over genes expressed
(max log2 CPM >= 1 on either side) in both members of the pair; the dot is the median over intestinal studies.
Intestinal studies: the Elmentaite 2020 + 2021 pair (fold changes averaged), Hickey 2023, Egozi 2021. The seven
Elmentaite 2020 child Crohn's donors are removed first (roster). A side with fewer than 20 cells is skipped.

Inputs (name=path):
  gut_sums gut_detected gut_groups   primary intestinal EEC pseudobulk by primary study label (POSITIVE = EC)
  sc_sums sc_detected sc_groups      SC-islet endocrine pseudobulk by primary study label
  ind_sums ind_detected ind_groups   induced EEC (Singh 2024) pseudobulk
  hio_sums hio_detected hio_groups   HIO graft EEC (McCauley 2023) pseudobulk, Endocrine_1 (EC) vs Endocrine_2
  roster                             intestinal reference roster (the excluded Elmentaite 2020 donors)
Outputs:
  pairwise       every population x intestinal study correlation (HIO graft excluded)
  dots           their median per population and study (the dots of Figure 3C)
  hio_pairwise   the same for the HIO graft (pooled and per organoid)
  hio_dots       their medians
  gut_folds      intestinal EC fold changes per intestinal study (read by gene_groups.py)

  python signature_correlations.py gut_sums=... ... dots=results/figure3/panel_c_dots.tsv ...
"""
import gzip
import json
import sys

import numpy as np
import pandas as pd

arg = dict(a.split("=", 1) for a in sys.argv[1:])
MIN_CELLS = 20
MAIN_STUDIES = ["GSE151117", "GSE167880", "GSE185036", "GSE199636", "GSE202497", "GSE226346"]

roster = json.load(open(arg["roster"], encoding="utf-8"))
E20, CROHNS_UNITS = roster["family_id"], tuple(roster["excluded_units"])
E21, HICKEY, EGOZI = "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC", "SF_EGOZI2021_DEVELOPMENTAL_EEC"
SOURCES = {"Elmentaite 2020 + 2021": [E20, E21], "Hickey 2023": [HICKEY], "Egozi 2021": [EGOZI]}


def load(tag):
    s = pd.read_csv(arg[tag + "_sums"], sep="\t", index_col="gene")
    d = pd.read_csv(arg[tag + "_detected"], sep="\t", index_col="gene")
    g = pd.read_csv(arg[tag + "_groups"], sep="\t", dtype={"unit": str})
    assert list(s.columns) == list(d.columns) == list(g.group), tag
    assert not s.index.has_duplicates, tag
    return s, d, g.set_index("group")


def pool(s, d, g, mask):
    cols = g.index[mask]
    return s[cols].sum(axis=1), d[cols].sum(axis=1), int(g.loc[cols, "n_cells"].sum())


def log2cpm(counts):
    return np.log2(counts / counts.sum() * 1e6 + 1)


def contrast(s, d, g, masks, a, b):
    sa, da, na = pool(s, d, g, masks[a]); sb, db, nb = pool(s, d, g, masks[b])
    if na < MIN_CELLS or nb < MIN_CELLS:
        return None
    la, lb = log2cpm(sa), log2cpm(sb)
    return {"lfc": la - lb, "la": la, "lb": lb, "fa": da / na, "fb": db / nb, "na": na, "nb": nb}


def average(cs):
    cs = [c for c in cs if c is not None]
    if not cs:
        return None
    genes = cs[0]["lfc"].index
    assert all(c["lfc"].index.equals(genes) for c in cs)
    out = {k: sum(c[k] for c in cs) / len(cs) for k in ("lfc", "la", "lb", "fa", "fb")}
    out.update(na=sum(c["na"] for c in cs), nb=sum(c["nb"] for c in cs), datasets=len(cs))
    return out


gut_s, gut_d, gut_g = load("gut")
sc_s, sc_d, sc_g = load("sc")
ind_s, ind_d, ind_g = load("ind")
hio_s, hio_d, hio_g = load("hio")

# roster: remove the Elmentaite 2020 Crohn's donor units
cols = gut_g.index[(gut_g.family_id == E20) & gut_g.unit.isin(CROHNS_UNITS)]
if set(gut_g.loc[cols, "unit"]) != set(CROHNS_UNITS):
    sys.exit("not every roster unit is in the intestinal pseudobulk")
gut_g = gut_g.drop(cols); gut_s = gut_s[gut_g.index]; gut_d = gut_d[gut_g.index]


def gut_contrast(fam):
    m = gut_g.family_id == fam
    return contrast(gut_s, gut_d, gut_g, {"EC": m & (gut_g["class"] == "POSITIVE"), "nonEC": m & (gut_g["class"] == "NEGATIVE")}, "EC", "nonEC")


fam_c = {f: gut_contrast(f) for f in [E20, E21, HICKEY, EGOZI]}
GUT = {name: average([fam_c[f] for f in fams]) for name, fams in SOURCES.items()}
if any(c is None for c in GUT.values()) or GUT["Elmentaite 2020 + 2021"]["datasets"] != 2:
    sys.exit("an intestinal study fell below the cell floor")

studies = sorted(sc_g.family_id.unique())
if studies != MAIN_STUDIES:
    sys.exit("the SC pseudobulk does not hold the six 10x studies")


def sc_contrast(study, pop):
    m = sc_g.family_id == study
    rest = (sc_g["class"] != pop) & ((sc_g["class"] != "SC_SIC") if pop != "SC_SIC" else True)
    return contrast(sc_s, sc_d, sc_g, {"A": m & (sc_g["class"] == pop), "B": m & rest}, "A", "B")


SC_POPS = ["SC_SIC", "SC_BETA", "SC_ALPHA"]
sc_c = {p: {st: c for st in studies if (c := sc_contrast(st, p)) is not None} for p in SC_POPS}
if sorted(sc_c["SC_SIC"]) != MAIN_STUDIES:
    sys.exit("an SC study has too few SC-SIC cells")

D6 = ind_g.unit == "D6"
IND = contrast(ind_s, ind_d, ind_g, {"A": D6 & (ind_g["class"] == "EC"), "B": D6 & (ind_g["class"] != "EC")}, "A", "B")

SPACE = sorted(set(gut_s.index) & set(sc_s.index))
SPACE_IDX = pd.Index(SPACE)


def expressed(c):
    return (np.maximum(c["la"], c["lb"]).reindex(SPACE_IDX) >= 1).fillna(False).astype(bool)


def pair(c1, c2):
    ok = expressed(c1) & expressed(c2)
    x, y = c1["lfc"].reindex(SPACE_IDX)[ok], c2["lfc"].reindex(SPACE_IDX)[ok]
    ok2 = x.notna() & y.notna()
    x, y = x[ok2], y[ok2]
    return float(x.rank().corr(y.rank())), int(len(x))


def expressed_half(cdict, genes):
    mx = pd.DataFrame({k: np.maximum(v["la"].reindex(genes), v["lb"].reindex(genes)) for k, v in cdict.items()})
    return (mx >= 1).sum(axis=1) >= np.ceil(mx.shape[1] / 2)


# intestinal EC fold changes per intestinal study (also used by the Figure 3F gene groups)
SIG_NAME = {"Elmentaite 2020 + 2021": "Gut Cell Atlas pair", "Hickey 2023": "Hickey 2023", "Egozi 2021": "Egozi 2021"}
GUT_LFC = pd.DataFrame({k: v["lfc"].reindex(SPACE) for k, v in GUT.items()})
GS = pd.DataFrame({"gut_ec_lfc": GUT_LFC.median(axis=1),
                   "gut_ec_up_in_half": (GUT_LFC >= 1).sum(axis=1) >= np.ceil(GUT_LFC.shape[1] / 2),
                   "gut_judgeable": expressed_half(GUT, SPACE)})
for src in GUT:
    GS[f"gut_lfc: {SIG_NAME[src]}"] = GUT_LFC[src]
GS["e20_alone_lfc"] = fam_c[E20]["lfc"].reindex(SPACE)
GS["e21_alone_lfc"] = fam_c[E21]["lfc"].reindex(SPACE)
GS.index.name = "gene"

if set(hio_g["class"]) != {"EC", "OTHER_EEC"}:
    sys.exit("HIO classes are not EC / OTHER_EEC")
HIO = {}
for name, um in (("McCauley 2023 both organoids", hio_g.unit.isin(["wt1", "wt2"])),
                 ("McCauley 2023 organoid wt1", hio_g.unit == "wt1"), ("McCauley 2023 organoid wt2", hio_g.unit == "wt2")):
    HIO[name] = contrast(hio_s, hio_d, hio_g, {"A": um & (hio_g["class"] == "EC"), "B": um & (hio_g["class"] == "OTHER_EEC")}, "A", "B")

POP_ROW = {"GUT_EC": "Gut EC vs other gut EEC", "SC_SIC": "SC-SIC vs other SC endocrine", "SC_BETA": "SC-beta vs other SC endocrine",
           "SC_ALPHA": "SC-alpha vs other SC endocrine", "INDUCED_EEC": "Induced EEC (EC vs other EEC)",
           "HIO_GRAFT_EC": "HIO graft EEC (EC vs other EEC)"}
members = [("GUT_EC", g, c) for g, c in GUT.items()]
for p in SC_POPS:
    members += [(p, st, c) for st, c in sc_c[p].items()]
members.append(("INDUCED_EEC", "Singh 2024 day 6", IND))
members += [("HIO_GRAFT_EC", k, c) for k, c in HIO.items()]
pw = []
for pop, study, c in members:
    for g, gc in GUT.items():
        if pop == "GUT_EC" and g == study:
            continue
        rho, n = pair(c, gc)
        pw.append({"population": pop, "study": study, "gut_study": g, "spearman": rho, "n_genes": n,
                   "n_cells_a": c["na"], "n_cells_b": c["nb"]})
PW = pd.DataFrame(pw)
DOTS = (PW.groupby(["population", "study"], sort=False)
        .agg(spearman_median=("spearman", "median"), n_pairs=("spearman", "size"), n_genes_min=("n_genes", "min"),
             n_genes_median=("n_genes", "median"), n_genes_max=("n_genes", "max")).reset_index())
DOTS["row"] = DOTS.population.map(POP_ROW)
hio = PW.population == "HIO_GRAFT_EC"
PW[~hio].reset_index(drop=True).to_csv(arg["pairwise"], sep="\t", index=False, lineterminator="\n")
DOTS[DOTS.population != "HIO_GRAFT_EC"].to_csv(arg["dots"], sep="\t", index=False, lineterminator="\n")
PW[hio].reset_index(drop=True).to_csv(arg["hio_pairwise"], sep="\t", index=False, lineterminator="\n")
DOTS[DOTS.population == "HIO_GRAFT_EC"].reset_index(drop=True).to_csv(arg["hio_dots"], sep="\t", index=False, lineterminator="\n")
with gzip.open(arg["gut_folds"], "wt", encoding="utf-8", newline="\n") as f:
    GS.to_csv(f, sep="\t", lineterminator="\n")
