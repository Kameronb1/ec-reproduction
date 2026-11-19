# Figure S6: tables behind the fetal pancreas panel. Scripts are in workflow/scripts/figure_s6/.
# The gene groups are Figure 2's result.

FS6_OK = str(RESULTS / "figure_s6_inputs_verified.txt")
FS6_S = script("figure_s6", "").rstrip("/\\")


rule figure_s6_scores:
    input:
        ok=FS6_OK,
        groups=result("figure2", "gene_groups.tsv.gz"),
        ol_counts=data("figure_s6/olaniru2023_pseudobulk_counts.tsv.gz"),
        ol_profiles=data("figure_s6/olaniru2023_profiles.tsv"),
        yu_counts=data("figure_s6/yu2021b_pseudobulk_counts.tsv.gz"),
        yu_profiles=data("figure_s6/yu2021b_profiles.tsv"),
    output:
        points=result("figure_s6", "points.tsv"),
        rows=result("figure_s6", "rows.tsv"),
    shell:
        '"{RSCRIPT}" "{FS6_S}/fetal_pancreas_scores.R" groups="{input.groups}" '
        'ol_counts="{input.ol_counts}" ol_profiles="{input.ol_profiles}" yu_counts="{input.yu_counts}" '
        'yu_profiles="{input.yu_profiles}" points="{output.points}" rows="{output.rows}"'
