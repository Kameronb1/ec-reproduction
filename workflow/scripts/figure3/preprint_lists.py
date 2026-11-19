"""Figure 3C: how many genes of each gene group rank near the top of the serotonin-producing cells' gene lists in two
published SC-islet studies.

  Liu 2025, Supplementary Table 5: 36 gene programs (rows) x 5,000 genes, z scores. Program h3 is the one the authors
      assign to SLC18A1+ TPH1+ cells. Rule: z > 2.
  Fiancette 2026, Table S3: marker genes per cluster for the H1 and H3 lines, no scores. Cluster 1 is the SLC18A1+
      cluster. Rule: among the first 500 listed (the authors used the first 500 of these lists for their GO analysis).
Each row also gives the best other program or cluster of the same study under the same rule, as a control.

Usage: preprint_lists.py groups=<gene groups .tsv.gz> liu=<csv> fiancette=<xlsx> out=<tsv>
"""
import csv, gzip, re, sys, zipfile
import xml.etree.ElementTree as ET

A = dict(a.split("=", 1) for a in sys.argv[1:])

with gzip.open(A["groups"], "rt", encoding="utf-8", newline="") as fh:
    groups = list(csv.DictReader(fh, delimiter="\t"))
GROUP = {"G3": "G3_SC_SIC_OWN_PROGRAM", "G2": "G2_EC_GENES_SC_SIC_HAVE", "G1": "G1_EC_GENES_SC_SIC_LACK"}
GS = {k: sorted({r["gene"] for r in groups if r["group"] == v}) for k, v in GROUP.items()}
assert [len(GS[k]) for k in GROUP] == [84, 34, 28], {k: len(v) for k, v in GS.items()}
NOTE = {"G2": "G2 includes TPH1, SLC18A1, FEV and LMX1A, markers the authors used to name these cells",
        "G1": "G1 genes were selected as low in SC-SICs; no intestinal cells in this dataset", "G3": ""}

# Liu 2025: z scores, program h3
with open(A["liu"], encoding="utf-8", newline="") as fh:
    rd = csv.reader(fh)
    genes = next(rd)[1:]
    Z = {r[0]: [float(v) for v in r[1:]] for r in rd}
assert len(Z) == 36 and len(genes) == 5000 and len(set(genes)) == 5000
# the 15 highest z scores of program h3 are the 15 genes the authors draw for h3 (their Supplementary Figure 7a)
FIG7A_H3 = {"BRINP3", "COL5A2", "AC007614.1", "HDAC9", "CBLN1", "ARPP21", "MME", "TPH1", "AL354771.1", "KCNJ3",
            "ADAMTSL3", "LMX1A", "KCNK13", "PIP5K1B", "ZMAT4"}
assert {genes[i] for i in sorted(range(5000), key=lambda i: -Z["h3"][i])[:15]} == FIG7A_H3
zpass = {p: {g for g, z in zip(genes, v) if z > 2} for p, v in Z.items()}


# Fiancette 2026: the workbook read with the standard library (sheet XML); cluster lists start on row 4
def read_xlsx(path):
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rid = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    z = zipfile.ZipFile(path)
    ss = ["".join(t.text or "" for t in si.iter("{%s}t" % ns["m"]))
          for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns)]
    rel = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    out = {}
    for s in ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", ns):
        target = rel[s.get(rid)].lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        rows = []
        for row in ET.fromstring(z.read(target)).iter("{%s}row" % ns["m"]):
            cells = {}
            for c in row.findall("m:c", ns):
                col = 0
                for ch in re.match(r"[A-Z]+", c.get("r")).group(0):
                    col = col * 26 + ord(ch) - 64
                v = c.find("m:v", ns)
                cells[col - 1] = (ss[int(v.text)] if c.get("t") == "s" else v.text) if v is not None else ""
            rows.append([(cells.get(i) or "").strip() for i in range(max(cells) + 1)] if cells else [])
        out[s.get("name")] = rows
    return out


wb = read_xlsx(A["fiancette"])
FSHEETS = [("H1", "TableS3b_SBC_H1"), ("H3", "TableS3c_SBC_H3")]
lists = {}
for line, sheet in FSHEETS:
    rows = wb[sheet]
    labels, clusters = rows[1], rows[2]
    assert labels[1].startswith("EC") and clusters[1] == "Cluster_1", (sheet, labels[:2], clusters[:2])
    for j, cl in enumerate(clusters):
        col = [r[j] for r in rows[3:] if j < len(r) and r[j]]
        if len(col) > 3:  # skips "cluster too small for marker"
            lists[(line, cl)] = (labels[j], col)

COLS = ["preprint", "cell_line", "gene_list", "rule", "gene_set", "genes_in_set", "genes_in_source_list",
        "genes_meeting_rule", "genes_meeting_rule_names", "best_other_list", "genes_meeting_rule_there", "note"]
rows = []
for k in GROUP:
    gs = set(GS[k])
    other = max((p for p in Z if p != "h3"), key=lambda p: (len(zpass[p] & gs), p))
    hit = sorted(zpass["h3"] & gs)
    rows.append(["Liu 2025 (preprint)", "H1 and HUES8 knockout village", "gene program h3 (SLC18A1+ TPH1+ cells)",
                 "z score > 2 (of 5,000 genes scored)", k, len(gs), len(gs & set(genes)), len(hit), ", ".join(hit),
                 "program " + other, len(zpass[other] & gs), NOTE[k]])
for line, _ in FSHEETS:
    for k in GROUP:
        gs = set(GS[k])
        lab, col = lists[(line, "Cluster_1")]
        first = set(col[:500])
        others = [(l, c) for (l, c) in lists if l == line and c != "Cluster_1"]
        best = max(others, key=lambda lc: (len(set(lists[lc][1][:500]) & gs), lc[1]))
        hit = sorted(first & gs)
        rows.append(["Fiancette 2026 (preprint)", line, "cluster 1 marker genes (SLC18A1+ cluster)",
                     "among the first 500 listed markers (lists carry no scores)", k, len(gs), len(set(col) & gs),
                     len(hit), ", ".join(hit), best[1].replace("_", " ").lower() + " (" + lists[best][0] + ")",
                     len(set(lists[best][1][:500]) & gs), NOTE[k]])

with open(A["out"], "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\t".join(COLS) + "\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
