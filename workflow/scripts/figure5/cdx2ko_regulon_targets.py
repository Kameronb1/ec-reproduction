"""Figure S17 (CDX2 regulon targets in the knockout): do CDX2's predicted regulon targets fall in Zhu 2023's CDX2 knockout
(GSE202493: one wild-type library MM844, one knockout library MM845)? Readouts, sign and null are those of Figure 5F
(cdx2ko_group_test.py); this script uses its own random-number stream, so its draws are independent of 5F's.

Sets (CDX2(+) regulons from network_regulons; reproducible = program_pass in all_programs):
  pooled: union of the reproducible CDX2(+) regulons' targets; "pooled outside G1-G3": the same without Figure 3 group genes;
  each CDX2(+) regulon alone, with and without G1-G3 genes.
Populations (transferred all-stage labels, the 12,615 nuclei of the study's GEO deposit): SC_EC (SC-SIC), SC_beta, SC_alpha,
and SC_EC minus SC_beta. Readouts, knockout minus wild type: pseudobulk log2 (equal depth by binomial thinning, 20 draws averaged) and
depth-stratified per-cell detection, as in 5F. Null: 500 random sets drawn gene by gene from the same expression bin (20
quantile bins of the population's pooled pseudobulk level plus a zero bin), excluding the set's genes and every named 5F set
gene (G1, G2, G3, beta and alpha identity). Primary test (fixed before any result): pooled set, SC_EC, pseudobulk, predicted
lower in the knockout (CDX2(+) activates its targets); one-sided share of random sets at or below. Everything else is
descriptive. Measurable target: detected (>= 1 UMI) in at least 5% of SC_EC nuclei (both genotypes pooled).
Targets named by an older symbol (the Hickey 2023 and Egozi 2021 networks; FAM49B is CYRIB in the knockout matrix) are
read under their current symbol from target_symbol_updates.tsv (43 of its 299 names occur in these regulons; each row's
Ensembl ID is checked against the matrix). Four names stay unmatched: HIF1A-AS2, LINC01481, MFSD14C, ZNF816-ZNF321P.

Context (descriptive): every regulon drawn in Figure 5E (program_values) scored in SC_EC with the pseudobulk readout and its
own 500 random sets, from a second random stream, and kept when it has at least 5 measurable targets. Regulons share
targets, so no test is made across them.

Seeds: 20261029 for the CDX2 sets (one numpy Generator, fixed call order), 20261030 for the context regulons.

Inputs:  counts=<recounted_counts.h5ad>  labels=<transferred_labels.tsv.gz>  groups=<gene_groups_definitions.tsv.gz>
         regulons=<network_regulons.tsv.gz>  all_programs=<all_programs.tsv>  program_values=<program_values.tsv>
         symbol_updates=<target_symbol_updates.tsv>
Outputs: effects=<set_effects.tsv>  subtracted=<beta_subtracted.tsv>  sets=<sets.tsv>  context=<context_regulons.tsv>
"""
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

a = dict(x.split("=", 1) for x in sys.argv[1:])

SEED, SEED_CTX, NRAND, REPS = 20261029, 20261030, 500, 20
POPS = ["SC_EC", "SC_beta", "SC_alpha"]
MEAS = 0.05                                  # measurable: detected in >= 5% of SC_EC nuclei
MIN_MEAS_CTX = 5
BETA_ID = ["INS", "IAPP", "PDX1", "NKX6-1", "MAFA"]
ALPHA_ID = ["GCG", "ARX", "IRX2", "TTR"]
PB, PC = "pseudobulk log2 (equal depth)", "per-cell detection (depth-stratified)"


def dec(x):
    return np.array([v.decode() if isinstance(v, bytes) else v for v in x])


def strat_diff(score, ko, lib):
    q = np.quantile(lib, np.linspace(0, 1, 11)); b = np.clip(np.searchsorted(q, lib, side="right") - 1, 0, 9)
    num, den = 0.0, 0
    for k in range(10):
        a_ = (b == k) & ko; c = (b == k) & ~ko
        if a_.sum() and c.sum():
            num += a_.sum() * (score[a_].mean() - score[c].mean()); den += a_.sum()
    return num / den if den else np.nan


with h5py.File(a["counts"], "r") as f:
    X = sparse.csr_matrix((f["X/data"][()], f["X/indices"][()], f["X/indptr"][()]), shape=tuple(f["X"].attrs["shape"]))
    cells = dec(f["obs/_index"][()]); genes = dec(f["var/_index"][()]); gene_ids = dec(f["var/gene_ids"][()])
assert len(set(genes)) == len(genes), "duplicate gene symbols"
T = pd.read_csv(a["labels"], sep="\t")
T = T[T["in_geo_deposit"].astype(str) == "True"].reset_index(drop=True)
pos = pd.Series(np.arange(len(cells)), index=cells)
X = X[pos[T["cell"]].values].tocsr()
gi = {g: i for i, g in enumerate(genes)}
gg = pd.read_csv(a["groups"], sep="\t")
G = {k: set(v) for k, v in gg.groupby("group")["gene"]}
GROUP_GENES = G["G1_EC_GENES_SC_SIC_LACK"] | G["G2_EC_GENES_SC_SIC_HAVE"] | G["G3_SC_SIC_OWN_PROGRAM"]
named = {gi[g] for g in GROUP_GENES | set(BETA_ID) | set(ALPHA_ID) if g in gi}

R = pd.read_csv(a["regulons"], sep="\t")
AP = pd.read_csv(a["all_programs"], sep="\t")
passed = dict(zip(AP["regulon"], AP["program_pass"].astype(str) == "True"))
U = pd.read_csv(a["symbol_updates"], sep="\t")
for old, new, eid in zip(U["network_symbol"], U["current_symbol"], U["ensembl_id"]):
    assert old not in gi and new in gi and gene_ids[gi[new]] == eid, (old, new, eid)
upd = dict(zip(U["network_symbol"], U["current_symbol"]))
tg = {r: [upd.get(g, g) for g in t.split(",") if g] for r, t in zip(R["regulon"], R["targets"].fillna(""))}
assert all(len(set(v)) == len(v) for v in tg.values()), "a regulon lists a target under both names"
renamed = {r: {upd[g] for g in t.split(",") if g in upd} for r, t in zip(R["regulon"], R["targets"].fillna(""))}
cdx = R[(R["tf"] == "CDX2") & (R["sign"] == "+")]["regulon"].tolist()
cdx = sorted(cdx, key=lambda r: (not passed[r], r))          # reproducible first
rep = [r for r in cdx if passed[r]]
assert len(rep) == 2, rep
pooled = sorted(set(tg[rep[0]]) | set(tg[rep[1]]))
setdefs = [("pooled reproducible", "pooled", "all targets", pooled),
           ("pooled reproducible", "pooled", "outside G1-G3", [g for g in pooled if g not in GROUP_GENES])]
for r in cdx:
    setdefs += [(("reproducible" if passed[r] else "not reproducible"), r, "all targets", tg[r]),
                (("reproducible" if passed[r] else "not reproducible"), r, "outside G1-G3", [g for g in tg[r] if g not in GROUP_GENES])]

ko_all = (T["genotype"] == "CDX2-KO").values
lib_all = np.asarray(X.sum(axis=1)).ravel()
cand_all = np.arange(len(genes))
rng = np.random.default_rng(SEED)
eff, bsub, per_pop = [], [], {}
for pop in POPS:
    idx = np.where((T["transferred_label_all_stages"] == pop).values)[0]
    Xp = X[idx]; ko = ko_all[idx]; lib = lib_all[idx]; det = (Xp > 0).tocsc()
    pk = np.asarray(Xp[ko].sum(axis=0)).ravel().astype(np.int64); pw = np.asarray(Xp[~ko].sum(axis=0)).ravel().astype(np.int64)
    n = min(pk.sum(), pw.sum())
    # REPS thinning draws, each gene's log2(1 + CPM) averaged over them
    lk, lw = np.zeros(len(pk)), np.zeros(len(pw))
    for _ in range(REPS):
        tk = rng.binomial(pk, n / pk.sum()) if pk.sum() > n else pk
        tw = rng.binomial(pw, n / pw.sum()) if pw.sum() > n else pw
        lk += np.log2(1 + tk / tk.sum() * 1e6) / REPS; lw += np.log2(1 + tw / tw.sum() * 1e6) / REPS
    pbd = lk - lw
    pooled_lvl = np.log2(1 + (pk + pw) / (pk + pw).sum() * 1e6)
    base = np.array([i for i in cand_all if i not in named])
    nz = base[pooled_lvl[base] > 0]; qs = np.quantile(pooled_lvl[nz], np.linspace(0, 1, 21))
    binof = lambda v: -1 if v <= 0 else int(np.clip(np.searchsorted(qs, v, side="right") - 1, 0, 19))
    bin_all = np.array([binof(v) for v in pooled_lvl])
    detfrac = np.asarray(det.mean(axis=0)).ravel()
    res = {}
    for kind, unit, scope, gl in setdefs:
        ii = np.array([gi[g] for g in gl if g in gi])
        excl = set(ii.tolist()) | named
        cand = np.array([i for i in base if i not in excl])
        pools = {b: cand[bin_all[cand] == b] for b in np.unique(bin_all[cand])}
        sc = np.asarray(det[:, ii].sum(axis=1)).ravel() / len(ii)
        obs_c, obs_p = strat_diff(sc, ko, lib), float(pbd[ii].mean())
        rsets = [np.array([rng.choice(pools[b]) for b in bin_all[ii]]) for _ in range(NRAND)]
        null_c = np.array([strat_diff(np.asarray(det[:, r].sum(axis=1)).ravel() / len(r), ko, lib) for r in rsets])
        null_p = np.array([pbd[r].mean() for r in rsets])
        key = (unit, scope)
        res[key] = {"ii": ii, "rsets": rsets, "obs_c": obs_c, "obs_p": obs_p, "null_c": null_c}
        for ro, obs, null in [(PC, obs_c, null_c), (PB, obs_p, null_p)]:
            eff.append({"population": pop, "set": unit, "scope": scope, "regulon_status": kind, "readout": ro,
                        "targets": len(gl), "in_matrix": len(ii),
                        "ko_minus_wt": obs, "random_median": float(np.median(null)), "random_q025": float(np.quantile(null, .025)),
                        "random_q975": float(np.quantile(null, .975)), "excess_over_random": obs - float(np.median(null)),
                        "random_share_at_or_below": float((np.sum(null <= obs) + 1) / (NRAND + 1)),
                        "random_share_at_or_above": float((np.sum(null >= obs) + 1) / (NRAND + 1))})
    per_pop[pop] = {"res": res, "pbd": pbd, "det": det, "ko": ko, "lib": lib, "detfrac": detfrac,
                    "bin_all": bin_all, "base": base}
    print(pop, "done", flush=True)

# measurable targets (SC_EC detection) and the share of them lower in the knockout, per population
E = per_pop["SC_EC"]
meas = E["detfrac"] >= MEAS
for row in eff:
    ii = per_pop[row["population"]]["res"][(row["set"], row["scope"])]["ii"]
    m = ii[meas[ii]]
    row["measurable_targets"] = int(len(m))
    row["share_measurable_lower_pseudobulk"] = float(np.mean(per_pop[row["population"]]["pbd"][m] < 0)) if len(m) else np.nan
    row["primary_test"] = (row["population"] == "SC_EC" and row["set"] == "pooled" and row["scope"] == "all targets"
                           and row["readout"] == PB)

B = per_pop["SC_beta"]
for (unit, scope), re_ in E["res"].items():
    ii = re_["ii"]
    obs = re_["obs_p"] - float(B["pbd"][ii].mean())
    null = np.array([E["pbd"][r].mean() - B["pbd"][r].mean() for r in re_["rsets"]])
    obs_c = re_["obs_c"] - B["res"][(unit, scope)]["obs_c"]
    null_c = np.array([re_["null_c"][k] - strat_diff(np.asarray(B["det"][:, r].sum(axis=1)).ravel() / len(r), B["ko"], B["lib"])
                       for k, r in enumerate(re_["rsets"])])
    for ro, o, nl in [(PC, obs_c, null_c), (PB, obs, null)]:
        bsub.append({"set": unit, "scope": scope, "readout": ro, "ec_minus_beta_change": o, "random_median": float(np.median(nl)),
                     "random_q025": float(np.quantile(nl, .025)), "random_q975": float(np.quantile(nl, .975)),
                     "excess_over_random": o - float(np.median(nl)),
                     "random_share_at_or_below": float((np.sum(nl <= o) + 1) / (NRAND + 1)),
                     "random_share_at_or_above": float((np.sum(nl >= o) + 1) / (NRAND + 1))})

# set membership, one row per target
rev = {v: k for k, v in upd.items()}
assert len(rev) == len(upd)
srows = []
for kind, unit, scope, gl in setdefs:
    srcs = rep if unit == "pooled" else [unit]
    for g in gl:
        j = gi.get(g)
        srows.append({"set": unit, "scope": scope, "gene": g,
                      "network_symbol": rev[g] if any(g in renamed[r] for r in srcs) else g, "in_matrix": j is not None,
                      "figure3_group": next((k.split("_")[0] for k, v in G.items() if g in v), ""),
                      "sc_ec_detected_share": float(E["detfrac"][j]) if j is not None else np.nan,
                      "measurable": bool(j is not None and meas[j]),
                      "sc_ec_pseudobulk_ko_minus_wt_log2": float(E["pbd"][j]) if j is not None else np.nan})

# context: every Figure 5E regulon, SC_EC pseudobulk only, its own random stream
PV = pd.read_csv(a["program_values"], sep="\t")
rng2 = np.random.default_rng(SEED_CTX)
ctx = []
for _, p in PV.iterrows():
    gl = tg[p["program"]]
    ii = np.array([gi[g] for g in gl if g in gi])
    nm = int(meas[ii].sum()) if len(ii) else 0
    row = {"factor": p["factor"], "regulon": p["program"], "network": p["network"], "reproducible": str(p["program_pass"]) == "True",
           "factor_direction": p["direction"], "targets": len(gl), "in_matrix": len(ii), "measurable_targets": nm,
           "scored": nm >= MIN_MEAS_CTX}
    if row["scored"]:
        excl = set(ii.tolist()) | named
        cand = np.array([i for i in E["base"] if i not in excl])
        pools = {b: cand[E["bin_all"][cand] == b] for b in np.unique(E["bin_all"][cand])}
        null = np.array([E["pbd"][np.array([rng2.choice(pools[b]) for b in E["bin_all"][ii]])].mean() for _ in range(NRAND)])
        obs = float(E["pbd"][ii].mean())
        row.update({"ko_minus_wt": obs, "random_median": float(np.median(null)), "excess_over_random": obs - float(np.median(null)),
                    "random_share_at_or_below": float((np.sum(null <= obs) + 1) / (NRAND + 1))})
    ctx.append(row)
C = pd.DataFrame(ctx)
sc = C["scored"]
C.loc[sc, "rank_excess_lowest_first"] = C.loc[sc, "excess_over_random"].rank(method="min")
C["scored_regulons"] = int(sc.sum())

pd.DataFrame(eff).to_csv(a["effects"], sep="\t", index=False, lineterminator="\n")
pd.DataFrame(bsub).to_csv(a["subtracted"], sep="\t", index=False, lineterminator="\n")
pd.DataFrame(srows).to_csv(a["sets"], sep="\t", index=False, lineterminator="\n")
C.to_csv(a["context"], sep="\t", index=False, lineterminator="\n")
print("context scored:", int(sc.sum()), "of", len(C), flush=True)
