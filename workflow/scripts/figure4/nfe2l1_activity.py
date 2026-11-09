"""Figure S16A, NFE2L1 bars: change in the NFE2L1 regulon target score (AUCell) in each graft comparison and each longer
in vitro culture comparison, for SC-SIC, SC-beta and SC-alpha.

The input holds, per comparison, population and arm (earlier = reference, later = test), the mean AUCell score of the
Augsornworawat 2020 NFE2L1 (+) regulon in each biological unit (mean over that unit's samples). Each arm's activity is
the mean over its units; the drawn value is test minus reference ("means" mode).

    python nfe2l1_activity.py units=<unit activity table> out=<table>
"""
import csv
import sys

import numpy as np

A = dict(a.split("=", 1) for a in sys.argv[1:])
REGULON = "GSE151117__NFE2L1(+)"
POPS = ["SC_SIC", "SC_BETA", "SC_ALPHA"]

with open(A["units"], newline="", encoding="utf-8") as handle:
    rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["regulon"] == REGULON and r["mode"] == "means"]
# the eight graft comparisons, then the six longer in vitro cultures, in the order of the S16A rows
comparisons = ["A-HUES8", "A-WS4", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-TXP",
               "B-w0-w3", "B-w0-w6", "Z-D21-D39", "V-w0-w5", "P-w2-m6", "P-w2-m12"]
assert set(comparisons) == {r["comparison"] for r in rows}

out = ["comparison\tpopulation\tregulon\tdifference"]
for c in comparisons:
    for p in POPS:
        arm = {}
        for which in ("reference", "test"):
            vals = [float(r["activity"]) for r in rows if r["comparison"] == c and r["population"] == p and r["arm"] == which]
            assert vals, (c, p, which)
            arm[which] = float(np.mean(vals))
        out.append("\t".join([c, p, REGULON, repr(arm["test"] - arm["reference"])]))   # shortest exact decimal
with open(A["out"], "w", encoding="utf-8", newline="\n") as handle:
    handle.write("\n".join(out) + "\n")
