"""Figure 5B (and the tables behind Figure S9A and S9B): candidate regulators with reproducible regulons.

What it computes. Each network (12 co-expression networks: six SC-islet, Elmentaite 2021 fetal and adult, Hickey 2023,
Egozi 2021 and two adult-islet networks) predicts regulons: a factor's predicted target set, signed (+) activated or (-)
repressed. Each regulon is scored (AUCell on its targets) in every dataset; the effect is the mean over samples of
"positive minus control" (SC-SIC minus SC-beta in SC-islet studies, EC minus other enteroendocrine cells in intestinal
studies, beta minus other endocrine in adult islets). This script applies the selection rule:
  1. SC-islets: SC studies the network was not trained on with >= 10 measured targets and >= 50% of the regulon's
     targets measured; at least 5 counted; the same nonzero sign in every counted study (direction d = sign of the
     median); matched two-sided p < 0.05 in every counted study.
  2. Intestine: three intestinal studies (Elmentaite 2020 + 2021 as one study, the equal-weight mean over its counted
     parts; Hickey 2023; Egozi 2021); every study the network was not trained on is counted and its effect has sign d.
  3. A regulon is reproducible (program_pass) when 1 and 2 hold. A factor (symbol and sign) is selected when its regulon
     is reproducible in >= 2 networks and in >= half of the networks that have it, all in one SC direction.
Adult islets are a check after selection, never an input. Also writes, for every regulon of every selected factor, the
drawn values (SC signed median, intestinal study means and their parts, adult mean).
Overlap of each regulon's targets with gene groups G1-G3 (expression-matched exact test, Benjamini-Hochberg) is an
annotation column. Targets named by an older symbol (the Hickey 2023 and Egozi 2021 networks; FAM49B is CYRIB) are
matched to the gene bins under their current symbol from target_symbol_updates.tsv; the selection reads the stored
AUCell results and does not use target names.

Inputs (name=path): regulons, training, effects, gene_bins, gene_groups, symbol_updates
Outputs: all_programs, factors, program_values, gut_parts
  python regulon_selection.py regulons=... training=... effects=... gene_bins=... gene_groups=... symbol_updates=... \
      all_programs=... factors=... program_values=... gut_parts=...
"""
import collections, csv, gzip, statistics, sys

a = dict(x.split("=", 1) for x in sys.argv[1:])

SC_STUDIES = ["GSE151117", "GSE167880", "GSE185036", "GSE199636", "GSE202497", "GSE226346"]
POOL = ["GSE151117", "GSE167880", "GSE185036", "GSE199636", "GSE202497", "GSE226346", "E-MTAB-9536",
        "support_small_adult_intestine", "HICKEY_2023_DRYAD", "EGOZI_2021_ZENODO", "GSE148073", "GSE221156"]
GUT_STUDIES = {  # study -> (parts, training datasets that remove the whole study)
    "ELMENTAITE": (["E-MTAB-9536__primary", "E-MTAB-9543__primary", "E-MTAB-8901__primary"], ["E-MTAB-9536", "E-MTAB-9543", "E-MTAB-8901"]),
    "HICKEY_2023": (["HICKEY_2023_DRYAD__primary"], ["HICKEY_2023_DRYAD"]),
    "EGOZI_2021": (["EGOZI_2021_ZENODO__primary"], ["EGOZI_2021_ZENODO"]),
}
ADULT = ["GSE148073", "GSE221156"]
MIN_SC, MIN_NETWORKS, MIN_FRACTION = 5, 2, 0.5


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


def sgn(x):
    return (x > 0) - (x < 0)


def write_tsv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def exact_matched(gs, targets, universe_bins):
    """Exact one-sided P(overlap >= observed) for a random gene set with gs's expression-bin composition."""
    from scipy.stats import hypergeom
    import numpy as np
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


# ---- inputs ----------------------------------------------------------------------------------------------------------
upd = {r["network_symbol"]: r["current_symbol"] for r in read(a["symbol_updates"])}
programs = {}
for r in read(a["regulons"]):
    assert r["regulon"] not in programs
    tg = [upd.get(g, g) for g in r["targets"].split(",")]
    assert len(set(tg)) == len(tg), "%s lists a target under both names" % r["regulon"]
    programs[r["regulon"]] = (r["network"], dict(tf=r["tf"], sign=r["sign"], targets=tg))
training = {r["network"]: r["training_studies"].split(",") for r in read(a["training"])}
eff = collections.defaultdict(dict)   # regulon -> target -> dict(effect, p, cov, external)
for r in read(a["effects"]):
    eff[r["regulon"]][r["target"]] = dict(effect=fnum(r["effect"]), p=fnum(r["matched_p"]), cov=r["coverage_pass"] == "TRUE",
                                          external=r["target"].replace("__primary", "") not in training[r["network"]])
sc = {reg: {t: x for t, x in v.items() if t in SC_STUDIES} for reg, v in eff.items()}

# ---- G-group overlap annotation ----------------------------------------------------------------------------------------
bins = {r["gene"]: r["expression_bin"] for r in read(a["gene_bins"])}
gg = read(a["gene_groups"])
G = {g: [r["gene"] for r in gg if r["group"] == g] for g in sorted({r["group"] for r in gg})}
assert all(x in bins for v in G.values() for x in v)
assert not set(upd) & set(bins), "an older symbol is itself in the gene bins"
ov = {}
for reg, (bank, v) in programs.items():
    t = set(v["targets"]) & set(bins)
    for grp, gs in G.items():
        k = sum(1 for g in gs if g in t)
        ov[(reg, grp)] = dict(overlap=k, targets=len(t), expected=len(gs) * len(t) / len(bins),
                              p=exact_matched(gs, t, bins) if t else None, genes=sorted(g for g in gs if g in t))
keys = [k for k in ov if ov[k]["p"] is not None]
for k, q in zip(keys, bh([ov[k]["p"] for k in keys])):
    ov[k]["q"] = q

# ---- regulon level -----------------------------------------------------------------------------------------------------
rows, detail = [], {}
for reg, (bank, v) in sorted(programs.items()):
    s = sc.get(reg, {})
    counted = {d: x for d, x in s.items() if x["cov"] and x["external"] and x["effect"] is not None}
    signs = {d: sgn(x["effect"]) for d, x in counted.items()}
    consistent = len(counted) >= MIN_SC and len(set(signs.values())) == 1 and 0 not in signs.values()
    d = (1 if statistics.median(x["effect"] for x in counted.values()) > 0 else -1) if counted else 0
    tr = set(training[bank])
    st_detail, eligible, gut_counted, gut_agree, study_means = [], 0, 0, 0, {}
    for st, (parts_, train_ds) in GUT_STUDIES.items():
        if tr & set(train_ds):
            st_detail.append(st + ":training"); continue
        eligible += 1
        parts = {t: eff.get(reg, {}).get(t) for t in parts_}
        cp = {t: x for t, x in parts.items() if x is not None and x["cov"] and x["effect"] is not None}
        if not cp:
            st_detail.append(st + ":not_counted"); continue
        gut_counted += 1
        study_mean = statistics.fmean(x["effect"] for x in cp.values())   # one averaged vote per study
        agree = sgn(study_mean) == d
        gut_agree += agree
        study_means[st] = study_mean
        st_detail.append("%s:%s(mean%+.5f;%s)" % (st, "agree" if agree else "disagree", study_mean,
                                                ",".join("%s%s" % (t.split("__")[0], "+" if sgn(x["effect"]) > 0 else "-" if sgn(x["effect"]) < 0 else "0")
                                                         for t, x in cp.items())))
    anchored = gut_agree == eligible   # every eligible intestinal study agrees
    ac = {t: x for t, x in eff.get(reg, {}).items() if t in ADULT and x["cov"] and x["external"] and x["effect"] is not None}
    detail[reg] = dict(gut_eligible=eligible, gut_detail=";".join(st_detail), **{"gut_mean__" + st: study_means.get(st) for st in GUT_STUDIES})
    beyond = sum(1 for x in counted.values() if x["p"] is not None and x["p"] < 0.05)
    beyond_all = bool(counted) and beyond == len(counted)
    row = dict(regulon=reg, bank=bank, tf=v["tf"], sign=v["sign"], sc_counted=len(counted),
               sc_direction={1: "SC-SIC higher", -1: "SC-SIC lower", 0: "NA"}[d],
               sc_same_sign=sum(1 for x in signs.values() if x == d), consistent=consistent,
               sc_beyond_random=beyond, sc_beyond_random_all=beyond_all,
               sc_median_abs_effect=statistics.median(abs(x["effect"]) for x in counted.values()) if counted else None,
               gut_counted=gut_counted, gut_agree=gut_agree, gut_anchored=anchored)
    row.update(program_pass=consistent and anchored and beyond_all,
               adult_scored=len(ac), adult_beta_like=sum(1 for x in ac.values() if sgn(x["effect"]) == -d))
    for grp in G:
        o = ov[(reg, grp)]; short = grp.split("_")[0]
        row.update({short + "_overlap": o["overlap"], short + "_expected": round(o["expected"], 3), short + "_q": o.get("q"),
                    short + "_genes": ",".join(o["genes"]) if o.get("q") is not None and o["q"] < 0.05 else ""})
    for dname in SC_STUDIES:
        x = s.get(dname)
        row["sc__" + dname] = None if x is None else x["effect"]
    row.update(detail[reg])
    rows.append(row)

# ---- factor level ------------------------------------------------------------------------------------------------------
def g23(r):
    return (r["G2_q"] if r["G2_q"] is not None else 1) < 0.05 or (r["G3_q"] if r["G3_q"] is not None else 1) < 0.05


assert max(collections.Counter((r["bank"], r["tf"], r["sign"]) for r in rows).values()) == 1
fac = collections.defaultdict(list)
for r in rows:
    fac[(r["tf"], r["sign"])].append(r)
factors = []
for (tf, sign), progs in fac.items():
    passing = [r for r in progs if r["program_pass"]]
    dirs = sorted({r["sc_direction"] for r in passing})
    selected = (len(passing) >= MIN_NETWORKS and len(passing) >= MIN_FRACTION * len(progs) and len(dirs) == 1)
    n_ov = sum(1 for r in passing if g23(r))
    ad = [r for r in passing if r["adult_scored"]]
    f = dict(
        factor=tf + "(" + sign + ")", tf=tf, sign=sign, networks_with_program=len(progs),
        networks_evaluable=sum(1 for r in progs if r["sc_counted"] >= MIN_SC), networks_passing=len(passing),
        passing_fraction=round(len(passing) / len(progs), 3), direction=dirs[0] if len(dirs) == 1 else ";".join(dirs),
        selected=selected,
        median_effect=statistics.median(r["sc_median_abs_effect"] for r in passing) if passing else None,
        passing_with_G2_or_G3_overlap=n_ov, passing_with_G2_overlap=sum(1 for r in passing if (r["G2_q"] or 1) < 0.05),
        passing_with_G3_overlap=sum(1 for r in passing if (r["G3_q"] or 1) < 0.05),
        show_in_5E=selected and n_ov > len(passing) / 2,
        adult_check_programs_scored=len(ad),
        adult_check_beta_like=sum(1 for r in ad if r["adult_beta_like"] == r["adult_scored"]))
    f.update(passing_networks=",".join(sorted(r["bank"] for r in passing)), networks=",".join(sorted(r["bank"] for r in progs)))
    factors.append(f)
factors.sort(key=lambda f: (not f["selected"], -f["networks_passing"], -(f["median_effect"] or 0), f["factor"]))
is_sel = {f["factor"]: f["selected"] for f in factors}
for r in rows:
    r["factor_selected"] = is_sel[r["tf"] + "(" + r["sign"] + ")"]

write_tsv(a["all_programs"], sorted(rows, key=lambda r: (not r["program_pass"], r["tf"], r["sign"], r["bank"])))
write_tsv(a["factors"], factors)

# ---- drawn values per regulon of every selected factor (read back from the written tables, as text) ---------------------
progs = {r["regulon"]: r for r in read(a["all_programs"])}
facs = [f for f in read(a["factors"]) if f["selected"] == "True"]
vrows, parts_rows = [], []
for f in facs:
    d_f = 1 if f["direction"] == "SC-SIC higher" else -1
    regs = [r for b in POOL for r, p in progs.items() if p["bank"] == b and p["tf"] + "(" + p["sign"] + ")" == f["factor"]]
    assert len(regs) == int(f["networks_with_program"])
    for reg in regs:
        p = progs[reg]; bank = p["bank"]
        counted = {k: x for k, x in sc.get(reg, {}).items() if x["cov"] and x["external"] and x["effect"] is not None}
        sc_med = statistics.median(x["effect"] for x in counted.values()) if counted else None
        ac = {t: x for t, x in eff.get(reg, {}).items() if t in ADULT and x["cov"] and x["external"] and x["effect"] is not None}
        ad_mean = statistics.fmean(x["effect"] for x in ac.values()) if ac else None
        tr = set(training[bank]); means = {}
        for stn, (parts_, train_ds) in GUT_STUDIES.items():
            if tr & set(train_ds):
                parts_rows.append(dict(program=reg, factor=f["factor"], network=bank, study=stn, part="", part_effect="",
                                       study_mean="", status="training"))
                continue
            cp = {t: eff.get(reg, {}).get(t) for t in parts_}
            cp = {t: x for t, x in cp.items() if x is not None and x["cov"] and x["effect"] is not None}
            m = statistics.fmean(x["effect"] for x in cp.values()) if cp else None
            means[stn] = m
            if not cp:
                parts_rows.append(dict(program=reg, factor=f["factor"], network=bank, study=stn, part="", part_effect="",
                                       study_mean="", status="not_counted"))
            for t, x in cp.items():
                parts_rows.append(dict(program=reg, factor=f["factor"], network=bank, study=stn, part=t.replace("__primary", ""),
                                       part_effect=x["effect"], study_mean=m, status="counted"))
        vrows.append(dict(factor=f["factor"], direction=f["direction"], d=d_f, program=reg, network=bank,
                          program_pass=p["program_pass"] == "True", own_dataset_is_adult=bank in ADULT,
                          sc_counted=len(counted), sc_signed_median=sc_med, sc_program_direction=p["sc_direction"],
                          sc_consistent=p["consistent"], sc_beyond_random_all=p["sc_beyond_random_all"],
                          gut_eligible=p["gut_eligible"], gut_counted=p["gut_counted"], gut_agree=p["gut_agree"],
                          gut_anchored=p["gut_anchored"],
                          **{"gut_mean__" + k: means.get(k) for k in GUT_STUDIES},
                          adult_datasets=";".join(sorted(ac)), adult_n=len(ac), adult_mean=ad_mean,
                          **{"adult__" + t: (ac[t]["effect"] if t in ac else None) for t in ADULT}))
write_tsv(a["program_values"], vrows)
write_tsv(a["gut_parts"], parts_rows)
print("selected factors:", len(facs), "regulons drawn:", len(vrows), "reproducible:", sum(r["program_pass"] for r in vrows))
