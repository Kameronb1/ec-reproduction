"""Shared helpers for agreement_per_donor.py and signature_scatter.py (Supplementary Figure S1B, S1C).

Loads the EC identity pseudobulk tables, removes the seven Elmentaite 2020 child donors with Crohn's disease from the
intestinal EC reference, and builds the pooled contrasts behind Figure 2A: per study and side, raw counts summed over
cells, log2(CPM + 1), population minus comparison (each side needs >= 20 cells); the Elmentaite 2020 + 2021 pair is one
intestinal source, the mean of its two datasets' fold changes.
"""
import numpy as np
import pandas as pd

E20, E21, HICKEY, EGOZI = ('SF_ELMENTAITE2020_DEVGUT', 'SF_ELMENTAITE2021_GUT_ATLAS', 'SF_HICKEY2023_ADULT_EEC',
                           'SF_EGOZI2021_DEVELOPMENTAL_EEC')
CROHNS_UNITS = ('T017', 'T019', 'T176', 'T189', 'T197', 'T202', 'T203')   # Elmentaite 2020 child donors with Crohn's disease
SOURCES = {'Elmentaite 2020 + 2021': [E20, E21], 'Hickey 2023': [HICKEY], 'Egozi 2021': [EGOZI]}
EXPECTED_SC = ['GSE151117', 'GSE167880', 'GSE185036', 'GSE199636', 'GSE202497', 'GSE226346']
SC_POPS = ['SC_SIC', 'SC_BETA', 'SC_ALPHA']
MIN_CELLS = 20


def load(pb, tag):
    s = pd.read_csv(f'{pb}/{tag}_sums.tsv.gz', sep='\t', index_col='gene')
    d = pd.read_csv(f'{pb}/{tag}_detected.tsv.gz', sep='\t', index_col='gene')
    g = pd.read_csv(f'{pb}/{tag}_groups.tsv', sep='\t', dtype={'unit': str})
    assert list(s.columns) == list(d.columns) == list(g.group), tag
    assert not s.index.has_duplicates, tag
    return s, d, g.set_index('group')


def pool(s, d, g, mask):
    cols = g.index[mask]
    return s[cols].sum(axis=1), d[cols].sum(axis=1), int(g.loc[cols, 'n_cells'].sum())


def log2cpm(counts):
    return np.log2(counts / counts.sum() * 1e6 + 1)


def contrast(s, d, g, masks, a, b):
    sa, da, na = pool(s, d, g, masks[a]); sb, db, nb = pool(s, d, g, masks[b])
    if na < MIN_CELLS or nb < MIN_CELLS:
        return None
    la, lb = log2cpm(sa), log2cpm(sb)
    return {'lfc': la - lb, 'la': la, 'lb': lb, 'fa': da / na, 'fb': db / nb, 'na': na, 'nb': nb}


def average(cs):
    cs = [c for c in cs if c is not None]
    if not cs: return None
    genes = cs[0]['lfc'].index
    assert all(c['lfc'].index.equals(genes) for c in cs)
    out = {k: sum(c[k] for c in cs) / len(cs) for k in ('lfc', 'la', 'lb', 'fa', 'fb')}
    out.update(na=sum(c['na'] for c in cs), nb=sum(c['nb'] for c in cs), datasets=len(cs))
    return out


def build(pb):
    """Everything Figure 2A's pooled comparison needs, on the intestinal reference without the Crohn's donors."""
    ns = {}
    gut_s, gut_d, gut_g = load(pb, 'gut_primary')
    sc_s, sc_d, sc_g = load(pb, 'sc')
    # the seven Crohn's donors make up 11 pseudobulk columns holding 10 EC and 15 other EEC cells
    before = gut_g.groupby(['family_id', 'class']).n_cells.sum()
    crohn_cols = gut_g.index[(gut_g.family_id == E20) & gut_g.unit.isin(CROHNS_UNITS)]
    removed = {k: int(v) for k, v in gut_g.loc[crohn_cols].groupby('class').n_cells.sum().items()}
    assert set(gut_g.loc[crohn_cols, 'unit']) == set(CROHNS_UNITS) and len(crohn_cols) == 11
    assert removed == {'POSITIVE': 10, 'NEGATIVE': 15}
    assert (int(before[(E20, 'POSITIVE')]), int(before[(E20, 'NEGATIVE')])) == (55, 140)
    keep = ~gut_g.index.isin(crohn_cols)
    gut_s, gut_d, gut_g = gut_s.loc[:, keep], gut_d.loc[:, keep], gut_g.loc[keep]
    after = gut_g.groupby(['family_id', 'class']).n_cells.sum()
    assert (int(after[(E20, 'POSITIVE')]), int(after[(E20, 'NEGATIVE')])) == (45, 125)
    assert after.drop(E20, level=0).equals(before.drop(E20, level=0))

    def gut_masks(fam, unit_pred=None):
        m = gut_g.family_id == fam
        if unit_pred is not None: m &= gut_g.unit.map(unit_pred).astype(bool)
        return {'EC': m & (gut_g['class'] == 'POSITIVE'), 'nonEC': m & (gut_g['class'] == 'NEGATIVE')}

    # EC / other EEC cells per study arm (Elmentaite 2020 child arm without the Crohn's donors)
    arms = {'E20 fetal': (E20, lambda u: u[0].isdigit(), (36, 107)), 'E20 child': (E20, lambda u: u.startswith('T'), (9, 18)),
            'E21 fetal': (E21, lambda u: u.startswith('F'), (276, 358)), 'E21 adult': (E21, lambda u: u.startswith('A'), (43, 56)),
            'Egozi fetal': (EGOZI, lambda u: u.endswith('w'), (119, 153)),
            'Egozi neonatal': (EGOZI, lambda u: u.endswith('_neo'), (51, 130)), 'Hickey': (HICKEY, lambda u: True, (780, 1149))}
    for arm, (fam, pred, want) in arms.items():
        m = gut_masks(fam, pred)
        got = (int(gut_g.loc[m['EC'], 'n_cells'].sum()), int(gut_g.loc[m['nonEC'], 'n_cells'].sum()))
        assert got == want, (arm, got, want)
    fam_c = {f: contrast(gut_s, gut_d, gut_g, gut_masks(f), 'EC', 'nonEC') for f in [E20, E21, HICKEY, EGOZI]}
    gut_c = {src: average([fam_c[f] for f in fams]) for src, fams in SOURCES.items()}
    assert all(c is not None for c in gut_c.values())

    SC_STUDIES = sorted(sc_g.family_id.unique()); assert SC_STUDIES == EXPECTED_SC, SC_STUDIES

    def sc_contrast(study, pop):
        m = sc_g.family_id == study
        rest = (sc_g['class'] != pop) & ((sc_g['class'] != 'SC_SIC') if pop != 'SC_SIC' else True)
        return contrast(sc_s, sc_d, sc_g, {'A': m & (sc_g['class'] == pop), 'B': m & rest}, 'A', 'B')

    sc_c = {p: {st: c for st in SC_STUDIES if (c := sc_contrast(st, p)) is not None} for p in SC_POPS}
    assert sorted(sc_c['SC_SIC']) == EXPECTED_SC

    hio_s, hio_d, hio_g = load(pb, 'hio_mccauley_endocrine1')
    assert set(hio_g['class']) == {'EC', 'OTHER_EEC'} and set(hio_g.unit) == {'wt1', 'wt2'} and int(hio_g.n_cells.sum()) == 598
    HIO = contrast(hio_s, hio_d, hio_g, {'A': hio_g['class'] == 'EC', 'B': hio_g['class'] == 'OTHER_EEC'}, 'A', 'B')
    assert (HIO['na'], HIO['nb']) == (322, 276), 'pooled HIO is 322 EC against 276 other EEC'

    SPACE = sorted(set(gut_s.index) & set(sc_s.index))
    ns.update(gut_s=gut_s, gut_g=gut_g, sc_s=sc_s, sc_g=sc_g, hio_s=hio_s, hio_g=hio_g,
              fam_c=fam_c, gut_c=gut_c, sc_c=sc_c, HIO=HIO, SPACE=SPACE)
    return ns
