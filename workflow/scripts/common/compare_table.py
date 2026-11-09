"""Compare a recreated table with the published one and write a JSON report.

Usage: compare_table.py <observed> <expected> <report.json>
                        [--mode exact|rows_unordered|numeric] [--keys a,b] [--rel-tolerance 0]
                        [--set-columns a,b]

Modes
  exact           identical header and identical cell text, row by row
  rows_unordered  identical header and identical rows, in any order
  numeric         identical header and non-numeric cells; numeric cells may differ by at most
                  --rel-tolerance relative to the expected value (0 means equal as numbers).
                  With --keys, rows are matched on those columns; otherwise by position.

Option for any mode
  --set-columns   comma-separated columns holding "/"-separated gene lists whose order carries no
                  meaning; these cells are compared as multisets (same genes, same counts).

Any tolerance other than 0 must be justified in expected/<figure>.tsv. The script exits
non-zero when the tables differ, so a mismatch stops the workflow rather than passing quietly.
"""
import argparse
import gzip
import json
import math
import sys


def read(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        lines = handle.read().replace("\r\n", "\n").rstrip("\n").split("\n")
    rows = [line.split("\t") for line in lines]
    return [cell.strip('"') for cell in rows[0]], [[cell.strip('"') for cell in r] for r in rows[1:]]


def as_number(text):
    if text in ("", "NA", "NaN", "nan"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("observed")
    parser.add_argument("expected")
    parser.add_argument("report")
    parser.add_argument("--mode", default="exact", choices=["exact", "rows_unordered", "numeric"])
    parser.add_argument("--keys", default="-")
    parser.add_argument("--rel-tolerance", type=float, default=0.0)
    parser.add_argument("--set-columns", default="-")
    args = parser.parse_args()

    obs_header, obs = read(args.observed)
    exp_header, exp = read(args.expected)
    set_columns = [] if args.set_columns == "-" else args.set_columns.split(",")
    for column in set_columns:
        for header, rows in ((obs_header, obs), (exp_header, exp)):
            if column in header:
                i = header.index(column)
                for r in rows:
                    r[i] = "/".join(sorted(r[i].split("/")))
    result = {"observed": args.observed, "expected": args.expected, "mode": args.mode, "set_columns": set_columns,
              "rel_tolerance": args.rel_tolerance, "rows_observed": len(obs), "rows_expected": len(exp)}
    differences = []

    if obs_header != exp_header:
        differences.append({"header": {"observed": obs_header, "expected": exp_header}})
    elif args.mode == "rows_unordered":
        if sorted(map(tuple, obs)) != sorted(map(tuple, exp)):
            differences.append({"rows": "different multiset of rows"})
    elif len(obs) != len(exp):
        differences.append({"row_count": [len(obs), len(exp)]})
    else:
        if args.keys != "-":
            idx = [exp_header.index(k) for k in args.keys.split(",")]
            key = lambda r: tuple(r[i] for i in idx)
            obs, exp = sorted(obs, key=key), sorted(exp, key=key)
        max_rel = 0.0
        for i, (o_row, e_row) in enumerate(zip(obs, exp)):
            for column, o, e in zip(exp_header, o_row, e_row):
                if o == e:
                    continue
                o_num, e_num = as_number(o), as_number(e)
                if args.mode == "numeric" and o_num is not None and e_num is not None:
                    if o_num == e_num:
                        rel = 0.0
                    elif e_num == 0:
                        rel = math.inf
                    else:
                        rel = abs(o_num - e_num) / abs(e_num)
                    max_rel = max(max_rel, rel)
                    if rel <= args.rel_tolerance:
                        continue
                if len(differences) < 20:
                    differences.append({"row": i + 1, "column": column, "observed": o, "expected": e})
                else:
                    result["more_differences"] = True
        result["max_relative_difference"] = max_rel

    result["status"] = "MATCH" if not differences else "DIFFERENT"
    result["differences"] = differences
    with open(args.report, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, indent=1)
    if differences:
        # Snakemake deletes the report of a failed job, so the first differences are printed here too.
        shown = "\n".join(f"  {d}" for d in differences[:5])
        sys.exit(f"{args.observed} differs from {args.expected}. First differences:\n{shown}")


if __name__ == "__main__":
    main()
