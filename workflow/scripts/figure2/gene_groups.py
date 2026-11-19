"""Figure 2D: three gene groups from the intestinal EC signature, the SC-SIC program and the Figure 2C fit.

  G1 EC genes low in SC-SICs:      up >= 2-fold in intestinal EC in all three intestinal studies, and higher in
                                   intestinal EC in 2C.
  G2 shared EC genes:              up in intestinal EC in all three intestinal studies, in the SC-SIC program, and not
                                   lower in SC-SICs in 2C.
  G3 separate SC-SIC program:      in the SC-SIC program, up in intestinal EC in at most one intestinal study, and
                                   higher in SC-SICs in 2C.
The SC-SIC program is the genes up >= 2-fold (log2 CPM, pooled counts) in SC-SICs against the other endocrine
cells with published labels in at least 5 of the 6 SC-islet studies. Genes up in exactly two intestinal studies
("mixed") are listed separately. Tested genes are those of the 2C fit. The group names in the output
(G1_EC_GENES_SC_SIC_LACK, G2_EC_GENES_SC_SIC_HAVE, G3_SC_SIC_OWN_PROGRAM) and the column names are table values.

Inputs (name=path):
  sc_sums sc_detected sc_groups     SC-islet endocrine pseudobulk by published label
  gut_sums gut_groups               primary intestinal EEC pseudobulk (to find intestinal studies with zero counts)
  roster                            intestinal reference roster (the excluded Elmentaite 2020 donors)
  gut_folds                         intestinal EC fold changes per intestinal study (signature_correlations.py)
  fit profiles                      2C levels table and pooled voom profiles (donor_fit.R)
Outputs:
  definitions   the G1-G3 genes with their values (Figure 2D)
  mixed         genes up in exactly two intestinal studies that would otherwise fall in G1 or G2
  universe      the tested genes (background for the Figure 2E enrichment)
"""
import gzip
import io
import json
import sys

import numpy as np
import pandas as pd
from scipy import special   # t.sf(x, df) is special.stdtr(df, -x)

arg = dict(a.split("=", 1) for a in sys.argv[1:])
MIN_CELLS = 20
K_SC = 5                                           # SC-SIC program: up in >= 5 of 6 SC studies
MAIN_STUDIES = ["GSE151117", "GSE167880", "GSE185036", "GSE199636", "GSE202497", "GSE226346"]
tru = lambda s: s.astype(str).isin(["True", "TRUE", "true"])


def load(tag):
    s = pd.read_csv(arg[tag + "_sums"], sep="\t", index_col="gene")
    d = pd.read_csv(arg[tag + "_detected"], sep="\t", index_col="gene")
    g = pd.read_csv(arg[tag + "_groups"], sep="\t", dtype={"unit": str}).set_index("group")
    assert list(s.columns) == list(d.columns) == list(g.index) and not s.index.has_duplicates, tag
    return s, d, g


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


# ---- SC signature fields and the SC-SIC program ----------------------------------------------------------------
sc_s, sc_d, sc_g = load("sc")
if sorted(sc_g.family_id.unique()) != MAIN_STUDIES:
    sys.exit("the SC pseudobulk does not hold the six 10x studies")
sc_c = {}
for pop in ("SC_SIC", "SC_BETA", "SC_ALPHA"):
    sc_c[pop] = {}
    for study in sorted(sc_g.family_id.unique()):
        m = sc_g.family_id == study
        rest = (sc_g["class"] != pop) & ((sc_g["class"] != "SC_SIC") if pop != "SC_SIC" else True)
        c = contrast(sc_s, sc_d, sc_g, {"A": m & (sc_g["class"] == pop), "B": m & rest}, "A", "B")
        if c is not None:
            sc_c[pop][study] = c

gf = pd.read_csv(arg["gut_folds"], sep="\t", index_col="gene")
SPACE = list(gf.index)                              # genes measured in both the intestinal and SC pseudobulks


def frame(cdict, key, genes):
    return pd.DataFrame({k: v[key].reindex(genes) for k, v in cdict.items()})


def expressed_half(cdict, genes):
    mx = pd.DataFrame({k: np.maximum(v["la"].reindex(genes), v["lb"].reindex(genes)) for k, v in cdict.items()})
    return (mx >= 1).sum(axis=1) >= np.ceil(mx.shape[1] / 2)


sig = pd.DataFrame(index=pd.Index(SPACE, name="gene"))
for pop in sc_c:
    sig[pop.lower() + "_lfc"] = frame(sc_c[pop], "lfc", SPACE).median(axis=1)
    sig[pop.lower() + "_judgeable"] = expressed_half(sc_c[pop], SPACE)
    sig[pop.lower() + "_studies_up1"] = (frame(sc_c[pop], "lfc", SPACE) >= 1).sum(axis=1)
# This table is written as text and read back with pandas' default float parser, which can move the
# last binary digit, so the values are identical to the published ones.
sig = pd.read_csv(io.StringIO(sig.to_csv(sep="\t")), sep="\t", index_col="gene")
SIC_ALL = pd.DataFrame({k: v["lfc"] for k, v in sc_c["SC_SIC"].items()})
sic_up = (SIC_ALL >= 1).sum(axis=1)
sic_up = sic_up[sic_up >= 4].sort_index()          # genes up in >= 4 of 6 studies; the program takes >= K_SC
PROGRAM = set(sic_up.index[sic_up >= K_SC])

# ---- intestinal status (three intestinal studies) ---------------------------------------------------------------
G1N, G2N, G3N = "G1_EC_GENES_SC_SIC_LACK", "G2_EC_GENES_SC_SIC_HAVE", "G3_SC_SIC_OWN_PROGRAM"
GUT3 = {"Elmentaite 2020 + 2021": "gut_lfc: Gut Cell Atlas pair", "Hickey 2023": "gut_lfc: Hickey 2023", "Egozi 2021": "gut_lfc: Egozi 2021"}
GUT_FAMILIES = {"Elmentaite 2020 + 2021": ["SF_ELMENTAITE2020_DEVGUT", "SF_ELMENTAITE2021_GUT_ATLAS"],
                "Hickey 2023": ["SF_HICKEY2023_ADULT_EEC"], "Egozi 2021": ["SF_EGOZI2021_DEVELOPMENTAL_EEC"]}
for c in ["gut_ec_lfc", "gut_ec_up_in_half", "gut_judgeable"] + list(GUT3.values()):
    sig[c] = gf[c]
gut = sig[list(GUT3.values())].set_axis(list(GUT3), axis=1)
gut_up = (gut >= 1).sum(axis=1)                    # an unmeasured study counts as not up
if set(gut_up.index[gut_up >= 2]) != set(sig.index[tru(sig.gut_ec_up_in_half)]):
    sys.exit("intestinal study counts do not reproduce gut_ec_up_in_half")
status = pd.Series(np.select([gut_up == 3, gut_up == 2], ["SHARED", "MIXED"], "NOT_RECURRENT"), index=gut_up.index)

roster = json.load(open(arg["roster"], encoding="utf-8"))
gs = pd.read_csv(arg["gut_sums"], sep="\t", index_col="gene")
gg = pd.read_csv(arg["gut_groups"], sep="\t", dtype={"unit": str}).set_index("group")
drop = (gg.family_id == roster["family_id"]) & gg.unit.isin(roster["excluded_units"])
gs, gg = gs.loc[:, (~drop).values], gg.loc[~drop]
zero = pd.DataFrame({s: gs[gg.index[gg.family_id.isin(f)]].sum(axis=1) == 0 for s, f in GUT_FAMILIES.items()})
zero_txt = zero.apply(lambda r: ";".join(c for c in zero.columns if r[c]), axis=1)

# SC-SIC study counts (studies with the gene up >= 2-fold): from the signature table for genes it can judge
# (expressed in at least half the studies); for genes up in >= 4 studies that are outside the signature table, from sic_up
sc_up1 = sig.sc_sic_studies_up1
judge = sig.sc_sic_judgeable.astype(bool)
SIC_COUNT = pd.concat([sc_up1.where(judge).astype("Int64"), sic_up[sic_up.index.difference(sc_up1.index)].astype("Int64")])


def sic_count(genes):
    return SIC_COUNT.reindex(pd.Index(genes)).array


# ---- the Figure 2C fit -----------------------------------------------------------------------------------------
fit = pd.read_csv(arg["fit"], sep="\t", low_memory=False).set_index("gene")
DF_STUDY = 7                                       # 3 intestinal + 6 SC study blocks - 2
p7 = 2 * special.stdtr(DF_STUDY, -fit.t.abs().values)
o = np.argsort(p7); q = p7[o] * len(p7) / np.arange(1, len(p7) + 1)
fdr7 = np.empty(len(p7)); fdr7[o] = np.minimum(np.minimum.accumulate(q[::-1])[::-1], 1)
fit["FDR_study_df7"] = pd.Series(fdr7, index=fit.index)
fit["robust_study_level"] = (fit.FDR_study_df7 < 0.05) & (fit.call != "NO_CALL")
TESTED = set(fit.index)

st = status.reindex(sorted(TESTED)).fillna("NOT_RECURRENT")
call = fit.call
G = {G1N: {g for g in TESTED if st[g] == "SHARED" and call[g] == "NATIVE_EC_HIGHER"},
     G2N: {g for g in TESTED & PROGRAM if st[g] == "SHARED" and call[g] in ("NO_CALL", "SC_SIC_HIGHER")},
     G3N: {g for g in TESTED & PROGRAM if st[g] == "NOT_RECURRENT" and call[g] == "SC_SIC_HIGHER"}}

prof = pd.read_csv(arg["profiles"], sep="\t").set_index("gene")
T = pd.DataFrame([{"group": grp, "gene": g} for grp in (G1N, G2N, G3N) for g in G[grp]])


def col(src, c, genes=None):
    return src[c].reindex(T.gene if genes is None else genes).values


T["panel5_logFC"] = col(fit, "logFC"); T["panel5_FDR"] = col(fit, "FDR")
T["panel5_core_welch"] = False
T["gut_ec_lfc"] = col(sig, "gut_ec_lfc"); T["sic_in_dish_lfc"] = col(sig, "sc_sic_lfc")
T["sic_in_dish_studies_up"] = sic_count(T.gene)
T["sc_beta_in_dish_lfc"] = col(sig, "sc_beta_lfc")
for s in prof.columns:
    T[f"level_{s}"] = col(prof, s)
T["paneld_call"] = col(fit, "call"); T["paneld_primary_level"] = col(fit, "primary_level"); T["paneld_sc_level"] = col(fit, "sc_level")
T["paneld_t"] = col(fit, "t"); T["paneld_FDR_study_df7"] = col(fit, "FDR_study_df7")
T["robust_study_level"] = col(fit, "robust_study_level").astype(bool)
for s in GUT3:
    T[f"gut_lfc_{s}"] = gut[s].reindex(T.gene).values
T["gut_up_of_3"] = gut_up.reindex(T.gene).fillna(0).astype(int).values
T["gut_status"] = st.reindex(T.gene).values
T["gut_zero_count_studies"] = zero_txt.reindex(T.gene).fillna("").values
T["sic_up_of_6"] = T.sic_in_dish_studies_up
T["g2_higher_in_sc_sic"] = (T.group == G2N) & (T.paneld_call == "SC_SIC_HIGHER")
T["g3_up_in_no_gut_study"] = (T.group == G3N) & (T.gut_up_of_3 == 0)
T["label_gene"] = tru(fit.label_gene.reindex(T.gene)).values
T = T.sort_values(["group", "paneld_t"], key=lambda s: s if s.name == "group" else -s.abs()).reset_index(drop=True)

mix = sorted(g for g in TESTED if st[g] == "MIXED")
M = pd.DataFrame({"gene": mix})
M["paneld_call"] = col(fit, "call", M.gene); M["in_program"] = M.gene.isin(PROGRAM)
M["would_be"] = np.select([M.paneld_call == "NATIVE_EC_HIGHER", M.in_program & M.paneld_call.isin(["NO_CALL", "SC_SIC_HIGHER"])], ["G1", "G2"], "")
M["would_be_g3_if_mixed_counted_not_shared"] = M.in_program & (M.paneld_call == "SC_SIC_HIGHER")
M = M[M.would_be != ""].copy()
for s in GUT3:
    M[f"gut_lfc_{s}"] = gut[s].reindex(M.gene).values
M["gut_study_not_up"] = [next(s for s in GUT3 if not gut.at[g, s] >= 1) for g in M.gene]
M["gut_zero_count_studies"] = zero_txt.reindex(M.gene).fillna("").values
M["paneld_logFC"] = col(fit, "logFC", M.gene); M["paneld_FDR"] = col(fit, "FDR", M.gene)
M["sic_up_of_6"] = sic_count(M.gene)
M = M.sort_values(["would_be", "gene"]).reset_index(drop=True)


def tsv_text(df):
    """LF, NA for missing, logicals as TRUE/FALSE."""
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == bool:
            df[c] = df[c].map({True: "TRUE", False: "FALSE"})
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False, lineterminator="\n", na_rep="NA")
    return buf.getvalue()


with gzip.open(arg["definitions"], "wt", encoding="utf-8", newline="\n") as f:
    f.write(tsv_text(T))
with open(arg["mixed"], "w", encoding="utf-8", newline="\n") as f:
    f.write(tsv_text(M))
with open(arg["universe"], "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(sorted(TESTED)) + "\n")
print({k: len(v) for k, v in G.items()}, "mixed", len(M))
