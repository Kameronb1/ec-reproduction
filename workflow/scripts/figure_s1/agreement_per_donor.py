"""Supplementary Figure S1B: does Figure 2A's agreement with intestinal EC hold when every donor or sample counts once?

For each population (intestinal EC, SC-SIC, SC-beta, SC-alpha, McCauley 2023 HIO graft EC) the fold change (population
vs comparison cells) is computed two ways and correlated (Spearman, over genes expressed at log2 CPM >= 1 on either side
in both members of a pair) with each intestinal EC source; each dot is the median over the intestinal sources.
  pooled    Figure 2A's pooled pseudobulk (all cells of a study and side summed; >= 20 cells per side).
  balanced  within each unit (intestine: donor; SC: sequenced sample; HIO: organoid)
            raw counts are summed per side, log2(CPM + 1), population minus comparison; a side needs >= 5 cells (10 and 20
            as sensitivity); the study's fold change is the mean over units that keep both sides, the Elmentaite pair the
            mean of its two datasets.
Intestinal EC reference: the seven Elmentaite 2020 child donors with Crohn's disease are excluded.

Inputs (name=path):
  lib       signature_lib.py
  pb        folder of the EC identity pseudobulk tables (<tag>_sums.tsv.gz, <tag>_detected.tsv.gz, <tag>_groups.tsv for
            gut_primary, sc, hio_mccauley_endocrine1)
  sc_counts folder of the SC endocrine cell counts (counts_indptr.i32, counts_indices.i32, counts_data.i32, genes.txt,
            metadata.tsv): genes x cells, compressed sparse columns, int32
Output:
  out       dots_pooled_vs_balanced.tsv (one row per dot)
"""
import importlib.util, sys
import numpy as np, pandas as pd, scipy.sparse as sp

A = dict(x.split("=", 1) for x in sys.argv[1:])
spec = importlib.util.spec_from_file_location('signature_lib', A['lib']); L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)
ns = L.build(A['pb'])
SPACE_IDX = pd.Index(ns['SPACE'])
log2cpm = L.log2cpm
MINS = (5, 10, 20)
MAIN_MIN = 5


def unit_balanced(sums, groups, unit_col, pos, neg, min_cells):
    """groups: DataFrame (index = column name in sums) with unit_col, 'class', 'n_cells'. Returns (study fc, unit table)."""
    fcs, rows = [], []
    for u, g in groups.groupby(unit_col, sort=True):
        a, b = g.index[g['class'] == pos], g.index[g['class'].isin(neg)]
        na, nb = int(g.loc[a, 'n_cells'].sum()), int(g.loc[b, 'n_cells'].sum())
        keep = na >= min_cells and nb >= min_cells
        rows.append({'unit': u, 'cells_population': na, 'cells_comparison': nb, 'kept': keep})
        if keep:
            fcs.append(log2cpm(sums[a].sum(axis=1)) - log2cpm(sums[b].sum(axis=1)))
    fc = pd.concat(fcs, axis=1).mean(axis=1) if fcs else None
    return fc, pd.DataFrame(rows)


# ---- SC per sequenced sample: sum each sample's cell counts per class ------------------------------------------------
SCP = A['sc_counts']
meta = pd.read_csv(f'{SCP}/metadata.tsv', sep='\t', dtype=str)
genes = open(f'{SCP}/genes.txt', encoding='utf-8').read().split('\n')
genes = [g for g in genes if g]
nC, nG = len(meta), len(genes)
indptr = np.fromfile(f'{SCP}/counts_indptr.i32', dtype='<i4')
assert len(indptr) == nC + 1
indices = np.fromfile(f'{SCP}/counts_indices.i32', dtype='<i4'); data = np.fromfile(f'{SCP}/counts_data.i32', dtype='<i4')
assert len(indices) == len(data) == indptr[-1]
X = sp.csc_matrix((data.astype(np.float64), indices, indptr), shape=(nG, nC)); del indices, data
key = meta.family_id + '|' + meta['sample'] + '|' + meta.target
codes, uniq = pd.factorize(key, sort=True)
ind = sp.csr_matrix((np.ones(nC), (np.arange(nC), codes)), shape=(nC, len(uniq)))
SC_SUMS = pd.DataFrame((X @ ind).toarray(), index=genes, columns=uniq); del X
parts = pd.Series(uniq).str.split('|', expand=True)
SC_GROUPS = pd.DataFrame({'family_id': parts[0].values, 'sample': parts[1].values, 'class': parts[2].values,
                          'n_cells': np.bincount(codes)}, index=uniq)
# the per-sample sums add up exactly to the pooled SC pseudobulk, study by study and class by class
sc_s, sc_g = ns['sc_s'], ns['sc_g']
assert list(sc_s.index) == genes
for (f, c), g in SC_GROUPS.groupby(['family_id', 'class']):
    ref = sc_s[sc_g.index[(sc_g.family_id == f) & (sc_g['class'] == c)]].sum(axis=1)
    assert np.array_equal(SC_SUMS[g.index].sum(axis=1).values, ref.values), (f, c)
    assert int(g.n_cells.sum()) == int(sc_g.loc[(sc_g.family_id == f) & (sc_g['class'] == c), 'n_cells'].sum()), (f, c)

gut_s, gut_g = ns['gut_s'], ns['gut_g']
hio_s, hio_g = ns['hio_s'], ns['hio_g']
GUT_FAMS = {'Elmentaite 2020 + 2021': [L.E20, L.E21], 'Hickey 2023': [L.HICKEY], 'Egozi 2021': [L.EGOZI]}
POPS = L.SC_POPS
SC_CLASSES = sorted(sc_g['class'].unique())
HIO_NAME = 'McCauley 2023 both organoids'


def build(min_cells):
    fc, units = {}, []
    for name, fams in GUT_FAMS.items():
        per = []
        for f in fams:
            g = gut_g[gut_g.family_id == f]
            x, u = unit_balanced(gut_s, g, 'unit', 'POSITIVE', ['NEGATIVE'], min_cells)
            units.append(u.assign(population='GUT_EC', study=name, dataset=f))
            if x is not None:
                per.append(x)
        fc[('GUT_EC', name)] = sum(per) / len(per) if per else None
    for pop in POPS:
        neg = [c for c in SC_CLASSES if c != pop and (pop == 'SC_SIC' or c != 'SC_SIC')]
        for f in sorted(SC_GROUPS.family_id.unique()):
            x, u = unit_balanced(SC_SUMS, SC_GROUPS[SC_GROUPS.family_id == f], 'sample', pop, neg, min_cells)
            units.append(u.assign(population=pop, study=f, dataset=f)); fc[(pop, f)] = x
    x, u = unit_balanced(hio_s, hio_g, 'unit', 'EC', ['OTHER_EEC'], min_cells)
    units.append(u.assign(population='HIO_GRAFT_EC', study=HIO_NAME, dataset='SF_GSE214852_HIO_GRAFT')); fc[('HIO_GRAFT_EC', HIO_NAME)] = x
    return fc, pd.concat(units, ignore_index=True)


GUT_POOLED = {name: L.average([ns['fam_c'][f] for f in fams]) for name, fams in GUT_FAMS.items()}
members = [('GUT_EC', g, c) for g, c in GUT_POOLED.items()]
for p in POPS:
    members += [(p, st, c) for st, c in ns['sc_c'][p].items()]
members.append(('HIO_GRAFT_EC', HIO_NAME, ns['HIO']))


def expressed(c):
    return (np.maximum(c['la'], c['lb']).reindex(SPACE_IDX) >= 1).fillna(False).astype(bool)


def rho(x1, x2, ok):
    if x1 is None or x2 is None:
        return np.nan, 0
    x, y = x1.reindex(SPACE_IDX)[ok], x2.reindex(SPACE_IDX)[ok]
    k = x.notna() & y.notna()
    return float(x[k].rank().corr(y[k].rank())), int(k.sum())


built = {m: build(m) for m in MINS}
rows = []
for pop, study, c in members:
    for g, gc in GUT_POOLED.items():
        if pop == 'GUT_EC' and g == study:
            continue
        ok = expressed(c) & expressed(gc)
        r = {'population': pop, 'study': study, 'gut_study': g}
        r['spearman_pooled'], r['n_genes'] = rho(c['lfc'], gc['lfc'], ok)
        for m in MINS:
            fc = built[m][0]
            r[f'spearman_balanced_min{m}'], _ = rho(fc[(pop, study)], fc[('GUT_EC', g)], ok)
        rows.append(r)
PW = pd.DataFrame(rows)
agg = {'pooled': ('spearman_pooled', 'median'), 'n_pairs': ('spearman_pooled', 'size'),
       'n_genes_min': ('n_genes', 'min'), 'n_genes_max': ('n_genes', 'max')}
agg.update({f'balanced_min{m}': (f'spearman_balanced_min{m}', 'median') for m in MINS})
DOTS = PW.groupby(['population', 'study'], sort=False).agg(**agg).reset_index()
for m in MINS:
    U = built[m][1]
    k = U.groupby(['population', 'study']).agg(units=('kept', 'size'), kept=('kept', 'sum')).reset_index()
    k[f'units_kept_min{m}'] = k.kept.astype(int).astype(str) + ' of ' + k.units.astype(str)
    DOTS = DOTS.merge(k[['population', 'study', f'units_kept_min{m}']], on=['population', 'study'], how='left')
kept5 = built[MAIN_MIN][1].groupby(['population', 'study']).kept.sum()
DOTS['flag_min5'] = [('NO_UNITS' if kept5[(p, s)] == 0 else 'ONE_UNIT' if kept5[(p, s)] == 1 else '') for p, s in zip(DOTS.population, DOTS.study)]
DOTS.to_csv(A['out'], sep='\t', index=False, lineterminator='\n')
