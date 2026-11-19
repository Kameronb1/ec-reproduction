"""Figure S9B: where the predicted targets of the selected candidate regulators (Figure 5B) fall.

What it computes. For each selected factor, its pooled targets = the union of the targets of its reproducible regulons
(one regulon per network).
  - Factors with targets higher in SC-SIC: within U = G2 + G3 (Figure 2's gene groups), the overlap of the pooled
    targets with G2, and with G3, against random groups of the same expression-bin composition drawn from U (exact
    expression-matched test, one-sided). Lean = toward G2 or toward G3 when that side's q < 0.05, else "no lean".
  - Factors with targets lower in SC-SIC: beta-side genes B = genes lower in SC-SIC than SC-beta (log2 fold
    change <= -1 of pooled-pseudobulk log2(CPM+1), >= 20 cells per side) in at least 5 SC-islet studies, minus the SC-SIC
    program (SC-SIC versus every other class, log2FC >= 1 in at least 5 studies). Overlap of the pooled targets with B by
    the exact matched test over the full gene universe, one-sided.
  - Factors whose pooled targets overlap with Jaccard >= 0.25 (same direction) are grouped into modules and tested too.
    Benjamini-Hochberg over every factor- and module-level test.
  - A chromatin annotation per factor with targets higher in SC-SIC from a G2/G3 motif screen (replicated family hits).
The tests of factors with targets lower in SC-SIC and of modules are written but not drawn.
Targets named by an older symbol (the Hickey 2023 and Egozi 2021 networks; FAM49B is CYRIB) are read under their current
symbol from target_symbol_updates.tsv before pooling.
Inputs (name=path): factors, all_programs, regulons, gene_bins, gene_groups, pseudobulk_sums, pseudobulk_groups,
                    motif_families, motif_hits, symbol_updates
Outputs: tests, factors_display
"""
import collections, csv, gzip, sys

import numpy as np
import pandas as pd

a = dict(x.split("=", 1) for x in sys.argv[1:])
JACCARD, MIN_CELLS = 0.25, 20
RECURRENCE_MIN = 5   # minimum number of SC-islet studies in which a gene's difference must recur


def read(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def exact_matched(gs, targets, universe_bins):
    """Exact one-sided P(overlap >= observed) for a random gene set with gs's expression-bin composition."""
    from scipy.stats import hypergeom
    by_bin = collections.defaultdict(list)
    for g, b in universe_bins.items():
        by_bin[b].append(g)
    prob = np.array([1.0])
    for k in sorted({universe_bins[g] for g in gs}):
        pop = by_bin[k]; draw = sum(universe_bins[g] == k for g in gs); hit = sum(1 for g in pop if g in targets)
        pp = hypergeom.pmf(np.arange(draw + 1), len(pop), hit, draw)
        prob = np.convolve(prob, pp)
    obs = sum(1 for g in gs if g in targets)
    return float(prob[obs:].sum())


def bh(ps):
    idx = sorted(range(len(ps)), key=lambda i: ps[i]); q = [None] * len(ps); m = len(ps); run = 1.0
    for rank in range(m - 1, -1, -1):
        i = idx[rank]; run = min(run, ps[i] * m / (rank + 1)); q[i] = run
    return q


def sc_contrasts(a_, b_):
    """Per SC study: log2(CPM+1) of pooled class sums, class a_ minus class b_ ('REST' = every other class); None when
    either side has < MIN_CELLS cells."""
    s = pd.read_csv(a["pseudobulk_sums"], sep="\t", index_col="gene")
    g = pd.read_csv(a["pseudobulk_groups"], sep="\t", dtype={"unit": str}).set_index("group")
    assert list(s.columns) == list(g.index)
    out = {}
    for st in sorted(g.family_id.unique()):
        m = g.family_id == st
        ma = m & (g["class"] == a_)
        mb = m & ((g["class"] != a_) if b_ == "REST" else (g["class"] == b_))
        na, nb = int(g.loc[ma, "n_cells"].sum()), int(g.loc[mb, "n_cells"].sum())
        if na < MIN_CELLS or nb < MIN_CELLS:
            continue
        sa, sb = s.loc[:, ma[ma].index].sum(axis=1), s.loc[:, mb[mb].index].sum(axis=1)
        out[st] = np.log2(sa / sa.sum() * 1e6 + 1) - np.log2(sb / sb.sum() * 1e6 + 1)
    return pd.DataFrame(out)


factors = read(a["factors"])
progs = read(a["all_programs"])
sel = [f for f in factors if f["selected"] == "True"]
upd = {r["network_symbol"]: r["current_symbol"] for r in read(a["symbol_updates"])}
banks = collections.defaultdict(dict)
for r in read(a["regulons"]):
    tg = [upd.get(g, g) for g in r["targets"].split(",")]
    assert len(set(tg)) == len(tg), "%s lists a target under both names" % r["regulon"]
    banks[r["network"]][r["regulon"]] = dict(tf=r["tf"], sign=r["sign"], targets=tg)

# pooled targets per selected factor (union over its replicated regulons)
pooled, passing_banks = {}, {}
for f in sel:
    bs = [r["bank"] for r in progs if r["tf"] == f["tf"] and r["sign"] == f["sign"] and r["program_pass"] == "True"]
    assert len(bs) == int(f["networks_passing"])
    t = set()
    for b in bs:
        v = [v for v in banks[b].values() if v["tf"] == f["tf"] and v["sign"] == f["sign"]]
        assert len(v) == 1; t |= set(v[0]["targets"])
    pooled[f["factor"]], passing_banks[f["factor"]] = t, bs
order = {f["factor"]: i for i, f in enumerate(sel)}
direction = {f["factor"]: f["direction"] for f in sel}

# modules: connected groups at Jaccard >= 0.25, within one direction
names = [f["factor"] for f in sel]
link = collections.defaultdict(set)
for i, x in enumerate(names):
    for y in names[i + 1:]:
        j = len(pooled[x] & pooled[y]) / len(pooled[x] | pooled[y])
        if direction[x] == direction[y] and j >= JACCARD:
            link[x].add(y); link[y].add(x)
module_of, modules = {}, []
for x in names:
    if x in module_of:
        continue
    comp, stack = set(), [x]
    while stack:
        z = stack.pop()
        if z not in comp:
            comp.add(z); stack.extend(link[z] - comp)
    members = sorted(comp, key=order.get)
    mid = "M%d" % (len(modules) + 1)
    for z in members:
        module_of[z] = mid
    modules.append(dict(module=mid, direction=direction[x], lead=members[0], members=members,
                        targets=set().union(*(pooled[z] for z in members))))

# gene universe and groups
bins = {r["gene"]: r["expression_bin"] for r in read(a["gene_bins"])}
gg = read(a["gene_groups"])
G2 = sorted({r["gene"] for r in gg if r["group"].startswith("G2_")})
G3 = sorted({r["gene"] for r in gg if r["group"].startswith("G3_")})
assert all(g in bins for g in G2 + G3) and not set(G2) & set(G3)
assert not set(upd) & set(bins), "an older symbol is itself in the gene bins"
U = set(G2) | set(G3)
binsU = {g: bins[g] for g in U}

# beta-side group
sic_rest = sc_contrasts("SC_SIC", "REST")
program = set(sic_rest.index[(sic_rest >= 1).sum(axis=1) >= RECURRENCE_MIN])
sic_beta = sc_contrasts("SC_SIC", "SC_BETA")
assert sic_beta.shape[1] >= RECURRENCE_MIN
B_all = set(sic_beta.index[(sic_beta <= -1).sum(axis=1) >= RECURRENCE_MIN]) - program
B = sorted(g for g in B_all if g in bins)

# chromatin families
fams = collections.defaultdict(set)
for r in read(a["motif_families"]):
    for m in (r["family_label"] + "/" + r["member_list"] + "/" + r["tf_name"]).replace("::", "/").split("/"):
        if m.strip():
            fams[r["family"]].add(m.strip().upper())
hits = read(a["motif_hits"])


def chromatin(tf):
    fs = {f for f, mem in fams.items() if tf.upper() in mem}
    if not fs:
        return "no chromatin test", ""
    rows = [h for h in hits if h["family"] in fs and h["group"] in ("G2", "G3")]
    rep = [h for h in rows if h["replicated_hit"] == "TRUE"]
    if rep:
        return "replicated motif enrichment", ";".join(sorted("%s %s (%s)" % (h["group"], h["window"], h["family_label"]) for h in rep))
    return "tested, not replicated", "min q " + ("%.3g" % min(float(h["max_q"]) for h in rows) if rows else "NA")


tests = []


def lean(label, T, level):
    t = (T & set(bins)) & U
    for grp, gs in (("G2", G2), ("G3", G3)):
        tests.append(dict(level=level, unit=label, test="lean_" + grp, targets_in_G2_G3=len(t),
                          overlap=len(t & set(gs)), plain_expected=round(len(t) * len(gs) / len(U), 2),
                          p=exact_matched(gs, t, binsU) if t else None))


def beta(label, T, level):
    t = T & set(bins)
    tests.append(dict(level=level, unit=label, test="beta_side", targets_in_G2_G3=None, overlap=len(t & set(B)),
                      plain_expected=round(len(t) * len(B) / len(bins), 2), p=exact_matched(B, t, bins) if t else None))


for f in names:
    (lean if direction[f] == "SC-SIC higher" else beta)(f, pooled[f], "factor")
for m in modules:
    if len(m["members"]) > 1:
        (lean if m["direction"] == "SC-SIC higher" else beta)(m["module"], m["targets"], "module")
idx = [i for i, t in enumerate(tests) if t["p"] is not None]
for i, q in zip(idx, bh([tests[i]["p"] for i in idx])):
    tests[i]["q"] = q
for t in tests:
    t.setdefault("q", None)
tq = {(t["unit"], t["test"]): t for t in tests}


def lean_call(u):
    q2, q3 = tq[(u, "lean_G2")]["q"], tq[(u, "lean_G3")]["q"]
    if q2 is not None and q2 < 0.05:
        return "toward G2 (EC genes SC-SIC have)"
    if q3 is not None and q3 < 0.05:
        return "toward G3 (SC-SIC's own part)"
    return "no lean"


rows = []
for f in sel:
    x = f["factor"]; up = direction[x] == "SC-SIC higher"; m = module_of[x]
    mrec = next(mm for mm in modules if mm["module"] == m)
    ch, chd = chromatin(f["tf"]) if up else ("no matching chromatin test", "")
    rows.append(dict(
        factor=x, direction=direction[x], module=m, module_members="+".join(mrec["members"]), module_lead=mrec["lead"],
        networks_passing=f["networks_passing"], networks_with_program=f["networks_with_program"],
        median_effect=f["median_effect"], pooled_targets=len(pooled[x]),
        g2_g3_lean=lean_call(x) if up else "",
        lean_G2_q=tq[(x, "lean_G2")]["q"] if up else None, lean_G3_q=tq[(x, "lean_G3")]["q"] if up else None,
        targets_in_G2_G3=tq[(x, "lean_G2")]["targets_in_G2_G3"] if up else None,
        beta_side_overlap=None if up else tq[(x, "beta_side")]["overlap"],
        beta_side_q=None if up else tq[(x, "beta_side")]["q"],
        chromatin=ch, chromatin_detail=chd, old_rule_D_show_in_5E=f["show_in_5E"],
        passing_networks=",".join(passing_banks[x])))

for name, data in (("factors_display", rows), ("tests", tests)):
    with open(a[name], "w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(data)
print("beta-side genes:", len(B_all), "in universe:", len(B), "tests:", len(tests))
