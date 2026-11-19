"""Figure 5C: what CDX2 knockout does to the gene groups of Figure 2 in SC-SIC, SC-beta and SC-alpha (Zhu 2023,
GSE202493: day-21 H1 SC-islets, one CDX2 wild-type library MM844 and one knockout library MM845).

Cell types are the labels transferred from Zhu 2023's reference with published labels (all-stage arm; see
cdx2ko_composition.py; SC_EC is SC-SIC). Genotype is the stored sample label (MM844 -> WT, MM845 -> CDX2-KO). Nuclei of one
genotype are pooled; the two groups are two libraries, so nothing here models biological replication: the random-gene
comparison tests whether a change is specific to a gene set, not whether it would replicate in another knockout
differentiation.

Counts are the study's reads recounted against GRCh38-2024-A (genes matched by its gene symbols, one per gene). Only the
12,615 nuclei that are also in the study's GEO deposit are used (in_geo_deposit in the label table).

Gene sets: G1, G2 and G3 from the gene-group table, with G2 split into the four EC markers (TPH1, SLC18A1, FEV, LMX1A;
set name "G2 serotonin core") and the other G2 genes; beta identity (INS, IAPP, PDX1, NKX6-1, MAFA); alpha identity
(GCG, ARX, IRX2, TTR). Readouts, knockout minus wild type within one transferred population (SC_EC, SC_beta, SC_alpha):
  per-cell detection: share of the set's genes detected (>= 1 UMI) per nucleus; difference taken within UMI deciles of
      that population, weighted by knockout nuclei per decile;
  pseudobulk: summed counts per genotype, the deeper side binomially thinned to the shallower total (20 draws, each
      gene's log2(1 + CPM) averaged over them), mean over the set of knockout minus wild type.
Each is compared with 500 random gene sets of the same size, drawn gene by gene from the same expression bin (20 quantile
bins of the population's pooled pseudobulk level plus a zero bin) of all 38,606 genes, excluding every named-set gene.
Every transferred nucleus is used (variant column: "all transferred"). SC_alpha is scored in a second loop, after SC_EC
and SC_beta, so their rows do not depend on it. Seed 20260930; one numpy Generator in a fixed call order.

Inputs:  counts=<recounted_counts.h5ad>  labels=<transferred_labels.tsv.gz>  groups=<gene-group table>
Outputs: effects=<group_effects.tsv>  nuclei=<nuclei.tsv>
"""
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

a = dict(x.split("=", 1) for x in sys.argv[1:])

SEED, NRAND, REPS = 20260930, 500, 20
POPS = ["SC_EC", "SC_beta"]
LATER_POPS = ["SC_alpha"]                    # second loop, after every SC_EC and SC_beta draw
SEROTONIN = ["TPH1", "SLC18A1", "FEV", "LMX1A"]
BETA_ID = ["INS", "IAPP", "PDX1", "NKX6-1", "MAFA"]
ALPHA_ID = ["GCG", "ARX", "IRX2", "TTR"]


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
    cells = dec(f["obs/_index"][()]); genes = dec(f["var/_index"][()])
assert len(set(genes)) == len(genes), "duplicate gene symbols"
T = pd.read_csv(a["labels"], sep="\t")
T = T[T["in_geo_deposit"].astype(str) == "True"].reset_index(drop=True)   # the 12,615 nuclei of the study's GEO deposit
pos = pd.Series(np.arange(len(cells)), index=cells)
X = X[pos[T["cell"]].values].tocsr()
gi = {}
for i, g in enumerate(genes):
    gi.setdefault(g, i)
gg = pd.read_csv(a["groups"], sep="\t")
G = {k: list(v) for k, v in gg.groupby("group")["gene"]}
sets0 = {"G1": G["G1_EC_GENES_SC_SIC_LACK"], "G2 serotonin core": SEROTONIN,
         "G2 other": [g for g in G["G2_EC_GENES_SC_SIC_HAVE"] if g not in SEROTONIN], "G3": G["G3_SC_SIC_OWN_PROGRAM"],
         "beta identity": BETA_ID, "alpha identity": ALPHA_ID}
missing = {k: [g for g in v if g not in gi] for k, v in sets0.items()}
print("genes not in the knockout matrix:", missing, flush=True)
sets = {k: np.array([gi[g] for g in v if g in gi]) for k, v in sets0.items()}
named = set(np.concatenate(list(sets.values())).tolist())
ko_all = (T["genotype"] == "CDX2-KO").values
lib_all = np.asarray(X.sum(axis=1)).ravel()
rng = np.random.default_rng(SEED)
eff, counts = [], []
passes = [("all transferred", np.ones(len(T), bool), POPS), ("all transferred", np.ones(len(T), bool), LATER_POPS)]
for vname, vmask, pops in passes:
    for pop in pops:
        idx = np.where(vmask & (T["transferred_label_all_stages"] == pop).values)[0]
        Xp = X[idx]; ko = ko_all[idx]; lib = lib_all[idx]; det = (Xp > 0).tocsc()
        counts.append({"variant": vname, "population": pop, "wt_nuclei": int((~ko).sum()), "ko_nuclei": int(ko.sum()),
                       "wt_median_umi": float(np.median(lib[~ko])), "ko_median_umi": float(np.median(lib[ko]))})
        pk = np.asarray(Xp[ko].sum(axis=0)).ravel().astype(np.int64); pw = np.asarray(Xp[~ko].sum(axis=0)).ravel().astype(np.int64)
        n = min(pk.sum(), pw.sum())
        # REPS thinning draws, each gene's log2(1 + CPM) averaged over them
        lk, lw = np.zeros(len(pk)), np.zeros(len(pw))
        for _ in range(REPS):
            tk = rng.binomial(pk, n / pk.sum()) if pk.sum() > n else pk
            tw = rng.binomial(pw, n / pw.sum()) if pw.sum() > n else pw
            lk += np.log2(1 + tk / tk.sum() * 1e6) / REPS; lw += np.log2(1 + tw / tw.sum() * 1e6) / REPS
        pbd = lk - lw
        pooled = np.log2(1 + (pk + pw) / (pk + pw).sum() * 1e6)
        cand = np.array([i for i in range(len(genes)) if i not in named and genes[i] == genes[gi[genes[i]]]])
        nz = cand[pooled[cand] > 0]; qs = np.quantile(pooled[nz], np.linspace(0, 1, 21))
        binof = lambda v: -1 if v <= 0 else int(np.clip(np.searchsorted(qs, v, side="right") - 1, 0, 19))
        cb = np.array([binof(v) for v in pooled[cand]]); pools = {b: cand[cb == b] for b in np.unique(cb)}
        for sname, ii in sets.items():
            sc = np.asarray(det[:, ii].sum(axis=1)).ravel() / len(ii)
            obs_c, obs_p = strat_diff(sc, ko, lib), float(pbd[ii].mean())
            bins = [binof(v) for v in pooled[ii]]
            rsets = [np.array([rng.choice(pools[b]) for b in bins]) for _ in range(NRAND)]
            null_c = np.array([strat_diff(np.asarray(det[:, r].sum(axis=1)).ravel() / len(r), ko, lib) for r in rsets])
            null_p = np.array([pbd[r].mean() for r in rsets])
            for ro, obs, null in [("per-cell detection (depth-stratified)", obs_c, null_c), ("pseudobulk log2 (equal depth)", obs_p, null_p)]:
                eff.append({"variant": vname, "population": pop, "set": sname, "readout": ro, "genes": len(ii),
                            "wt_score_mean": float(sc[~ko].mean()) if ro.startswith("per-cell") else float(lw[ii].mean()),
                            "ko_score_mean": float(sc[ko].mean()) if ro.startswith("per-cell") else float(lk[ii].mean()),
                            "ko_minus_wt": obs, "random_median": float(np.median(null)), "random_q025": float(np.quantile(null, 0.025)),
                            "random_q975": float(np.quantile(null, 0.975)), "excess_over_random": obs - float(np.median(null)),
                            "random_share_at_or_below": float((np.sum(null <= obs) + 1) / (NRAND + 1)),
                            "random_share_at_or_above": float((np.sum(null >= obs) + 1) / (NRAND + 1))})
    print(vname, " ".join(pops), "done", flush=True)
pd.DataFrame(eff).to_csv(a["effects"], sep="\t", index=False)
pd.DataFrame(counts).to_csv(a["nuclei"], sep="\t", index=False)
