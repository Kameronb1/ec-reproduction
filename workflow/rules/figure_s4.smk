# Figure S4: tables behind panels A and B. Scripts are in workflow/scripts/figure_s4/.
# The gene groups are Figure 2's result; the panels use Figure 3B's color scale (Figure 3's table g1_pooled_levels.tsv).

FS4_OK = str(RESULTS / "figure_s4_inputs_verified.txt")
FS4_S = script("figure_s4", "").rstrip("/\\")
FS4_PB = "shared/ec_identity_pseudobulk/"


rule figure_s4_pooled_levels:
    input:
        ok=FS4_OK,
        groups=result("figure2", "gene_groups.tsv.gz"),
        gut_sums=data(FS4_PB + "gut_primary_sums.tsv.gz"),
        gut_groups=data(FS4_PB + "gut_primary_groups.tsv"),
        hio_sums=data(FS4_PB + "hio_mccauley_endocrine1_sums.tsv.gz"),
        hio_groups=data(FS4_PB + "hio_mccauley_endocrine1_groups.tsv"),
        sc_sums=data(FS4_PB + "sc_sums.tsv.gz"),
        sc_groups=data(FS4_PB + "sc_groups.tsv"),
        gut_roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
    output:
        levels=result("figure_s4", "g2_g3_pooled_levels.tsv"),
        counts=result("figure_s4", "column_counts.tsv"),
        order=result("figure_s4", "gene_order.tsv"),
    shell:
        '"{RSCRIPT}" "{FS4_S}/pooled_levels.R" groups="{input.groups}" gut_sums="{input.gut_sums}" '
        'gut_groups="{input.gut_groups}" hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" sc_sums="{input.sc_sums}" '
        'sc_groups="{input.sc_groups}" gut_roster="{input.gut_roster}" levels="{output.levels}" counts="{output.counts}" '
        'order="{output.order}"'
