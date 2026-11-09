"""Figure S7B: gene by gene, how closely do each population's fold changes follow the intestinal EC signature?

Intestinal EC signature per gene (read from gut_source_fold_changes): median of the three intestinal sources'
EC-vs-other-EEC log2 fold changes (Figure 3C's pooled contrasts; Elmentaite 2020 + 2021 one source; the seven
Elmentaite 2020 child donors with Crohn's disease excluded).
A gene is judgeable in the intestine when log2 CPM >= 1 on either side in at least half the sources, and an intestinal EC
signature gene when it is up >= 2-fold in at least 2 of the 3 sources. Each of SC-SIC, SC-beta, SC-alpha, day-6 induced EEC (their
signatures from q1_gene_signatures: median over studies) and McCauley 2023 HIO graft EC (pooled, 322 EC vs 276 other
EEC) is compared with it: Spearman over genes judgeable on both sides, with and without the 18 marker genes; slope;
share of intestinal EC signature genes up >= 1.4-fold (log2 >= 0.5) in the population.

Inputs (name=path):
  lib        s3_lib.py
  pb         folder of the EC identity pseudobulk tables (only hio_mccauley_endocrine1_* is read here)
  gut        gut_source_fold_changes.tsv.gz (per gene: intestinal EC signature, up-in-half and judgeable flags)
  signatures q1_gene_signatures.tsv.gz (per-gene signatures of SC-SIC, SC-beta, SC-alpha, induced EEC; marker flags)
Outputs:
  out_genes  gene_signatures.tsv.gz (per gene: intestinal EC and the five populations' fold changes and judgeable flags)
  out_concordance  concordance.tsv
"""
import importlib.util, sys
import numpy as np, pandas as pd

A = dict(x.split("=", 1) for x in sys.argv[1:])
spec = importlib.util.spec_from_file_location('s3_lib', A['lib']); L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)

G = pd.read_csv(A['signatures'], sep='\t', index_col=0)
for c in [c for c in G.columns if c.endswith(('_judgeable', '_up_in_half')) or c == 'excluded18']:
    G[c] = G[c].astype(str).eq('True')
assert int(G.excluded18.sum()) == 18, 'q1_gene_signatures carries the 18 marker genes'
SPACE = G.index

GN = pd.read_csv(A['gut'], sep='\t', index_col=0)
assert GN.index.equals(SPACE), 'gut_source_fold_changes and q1_gene_signatures list the same genes'
GUT = GN[['gut_ec_lfc', 'gut_ec_up_in_half', 'gut_judgeable']].copy()
for c in ('gut_ec_up_in_half', 'gut_judgeable'): GUT[c] = GUT[c].astype(str).eq('True')
GUT_SIG, GUT_UP, GUT_OK = GUT.gut_ec_lfc, GUT.gut_ec_up_in_half, GUT.gut_judgeable
NOT_EXCL = ~G.excluded18


def concord(sig, ok_other, genes_mask=None, gut_sig=GUT_SIG, gut_ok=GUT_OK, gut_up=GUT_UP):
    ok = gut_ok & ok_other & sig.notna() & gut_sig.notna()
    if genes_mask is not None: ok &= genes_mask
    x, y = gut_sig[ok], sig[ok]
    rho = x.rank().corr(y.rank())
    slope = np.polyfit(x.to_numpy(), y.to_numpy(), 1)[0]
    up = gut_up & ok
    direction = float((sig[up] >= 0.5).mean()) if up.any() else np.nan
    return {'spearman': rho, 'slope': slope, 'direction_frac': direction, 'n_genes': int(ok.sum()), 'n_gut_up': int(up.sum())}


POPS = {'SC_SIC': 'sc_sic', 'SC_BETA': 'sc_beta', 'SC_ALPHA': 'sc_alpha', 'INDUCED_EEC': 'induced_eec'}
SIGS = {p: (G[f'{k}_lfc'], G[f'{k}_judgeable']) for p, k in POPS.items()}
hio_s, hio_d, hio_g = L.load(A['pb'], 'hio_mccauley_endocrine1')
assert set(hio_g['class']) == {'EC', 'OTHER_EEC'} and int(hio_g.n_cells.sum()) == 598
H = L.contrast(hio_s, hio_d, hio_g, {'A': hio_g['class'] == 'EC', 'B': hio_g['class'] == 'OTHER_EEC'}, 'A', 'B')
assert (H['na'], H['nb']) == (322, 276), 'pooled HIO is 322 EC against 276 other EEC'
SIGS['HIO_GRAFT_EC'] = (H['lfc'].reindex(SPACE), (np.maximum(H['la'], H['lb']).reindex(SPACE) >= 1).fillna(False).astype(bool))

rows = []
for p, (s, ok) in SIGS.items():
    for genes, mask in [('with the 18 excluded genes', None), ('without them', NOT_EXCL)]:
        r = concord(s, ok, mask); r.update(population=p, genes=genes); rows.append(r)
CON = pd.DataFrame(rows)[['population', 'genes', 'spearman', 'slope', 'direction_frac', 'n_genes', 'n_gut_up']]
CON.to_csv(A['out_concordance'], sep='\t', index=False, lineterminator='\n')
genes = pd.DataFrame({'gut_ec_lfc': GUT_SIG, 'gut_ec_up_in_half': GUT_UP, 'gut_judgeable': GUT_OK, 'excluded18': G.excluded18})
for p, (s, ok) in SIGS.items(): genes[f'{p}_lfc'] = s; genes[f'{p}_judgeable'] = ok
genes.to_csv(A['out_genes'], sep='\t', lineterminator='\n', compression={'method': 'gzip', 'mtime': 0})
