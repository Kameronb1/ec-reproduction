"""Figure S7B (shares): of the intestinal EC signature genes (up >= 2-fold in at least 2 of the 3 intestinal sources),
the share that is also up >= 2-fold (log2 >= 1) in each population; the 1.4-fold (log2 >= 0.5) share is written
alongside. Two gene sets, both without the 18 marker genes: each population's own judgeable genes ('own'), and as a check
the genes judgeable in intestinal EC and in all five populations at once ('common').

Inputs (name=path):
  genes  gene_signatures.tsv.gz written by s3_06_signature_scatter_five.py
Output:
  out    share_two_fold.tsv
"""
import sys
import pandas as pd

A = dict(x.split("=", 1) for x in sys.argv[1:])
G = pd.read_csv(A['genes'], sep='\t', index_col=0)
for c in [c for c in G.columns if c.endswith(('_judgeable', '_up_in_half')) or c == 'excluded18']: G[c] = G[c].astype(str).eq('True')
POPS = ['SC_SIC', 'SC_BETA', 'SC_ALPHA', 'INDUCED_EEC', 'HIO_GRAFT_EC']
base = G.gut_judgeable & G.gut_ec_lfc.notna() & ~G.excluded18
COMMON = base.copy()
for p in POPS: COMMON &= G[f'{p}_judgeable'] & G[f'{p}_lfc'].notna()
rows = []
for p in POPS:
    ok = base & G[f'{p}_judgeable'] & G[f'{p}_lfc'].notna()
    for gs, mask in (('own', ok), ('common', ok & COMMON)):
        up = G.gut_ec_up_in_half & mask
        y = G.loc[up, f'{p}_lfc']
        share_1p4, share_2 = float((y >= 0.5).mean()), float((y >= 1).mean())
        rows.append({'population': p, 'gene_set': gs, 'n_gut_ec_up': int(up.sum()), 'n_also_up_2fold': int((y >= 1).sum()),
                     'share_2fold': share_2, 'share_1p4fold': share_1p4})
pd.DataFrame(rows).to_csv(A['out'], sep='\t', index=False, lineterminator='\n')
