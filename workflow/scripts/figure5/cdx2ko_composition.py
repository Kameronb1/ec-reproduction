"""Figure S17F: cell-type shares in the CDX2 wild-type and knockout SC-islets (Zhu 2023, GSE202493).

The label table covers every recounted nucleus passing QC; only the 12,615 nuclei that are also in the study's GEO
deposit are counted here (in_geo_deposit). Each carries a cell-type label transferred from Zhu 2023's own SC-islet
reference (GSE202497, primary study labels). This script
  1. counts nuclei per genotype and transferred label, for both reference arms (all stages, day 21), exactly as the
     label transfer assigned them (genotype_by_label.tsv);
  2. compares the all-stage and day-21 shares of the five cell types Zhu's Figure 5D draws (ENP3, SC-alpha, SC-EC,
     SC-beta, SC-delta; shares rescaled to sum to 1 within each genotype) with the published bar heights, read from the
     paper's Figure 5D bar rectangles (figure5d_comparison.tsv). Genotype is the sample label (MM844 = WT, MM845 =
     CDX2-KO): one library per genotype, so the shares describe two samples, not replicated groups.

Inputs:  labels=<transferred_labels.tsv.gz>  bars=<published Figure 5D bar segments, JSON>
Outputs: by_label=<genotype_by_label.tsv>  comparison=<figure5d_comparison.tsv>
"""
import json
import sys

import pandas as pd

a = dict(x.split("=", 1) for x in sys.argv[1:])

qc = pd.read_csv(a["labels"], sep="\t")
qc = qc[qc["in_geo_deposit"].astype(str) == "True"].reset_index(drop=True)   # the 12,615 nuclei of the study's GEO deposit
arms = ["all_stages", "d21"]

# 1. genotype by transferred label
rows_g = []
for arm in arms:
    col = "transferred_label_" + arm
    for gt, d in qc.groupby("genotype"):
        vc = d[col].value_counts()
        for l, n in vc.items():
            rows_g.append({"arm": arm, "genotype": gt, "label": l, "nuclei": int(n), "fraction_of_genotype": round(n / len(d), 4),
                           "median_confidence": round(float(d.loc[d[col] == l, "transferred_confidence_" + arm].median()), 3)})
pd.DataFrame(rows_g).to_csv(a["by_label"], sep="\t", index=False)

# 2. the five drawn types against the published bars (bar height in points / axis height 65.36 pt)
FIVE = ["ENP3", "SC_alpha", "SC_EC", "SC_beta", "SC_delta"]
AXIS_PT = 313.57 - 248.21
bars = json.load(open(a["bars"], encoding="utf-8"))["segments_pt"]
published = {gt: {lab: (hi - lo) / AXIS_PT for lab, lo, hi in seg} for gt, seg in bars.items()}
rows = []
for arm in arms:
    col = "transferred_label_" + arm
    for gt in ["WT", "CDX2-KO"]:
        n = qc.loc[qc["genotype"] == gt, col].value_counts()
        tot = sum(int(n.get(l, 0)) for l in FIVE)
        for l in FIVE:
            ours = round(int(n.get(l, 0)) / tot, 3)
            zhu = round(published[gt][l], 3)
            rows.append({"arm": arm, "genotype": gt, "label": l, "ours": ours, "zhu_fig5d": zhu, "diff": round(ours - zhu, 3)})
pd.DataFrame(rows).to_csv(a["comparison"], sep="\t", index=False)
