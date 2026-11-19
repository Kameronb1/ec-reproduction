# Figure 3: tables behind panels A-C. Scripts are in workflow/scripts/figure3/.
# The gene groups are Figure 2's result.

F3_OK = str(RESULTS / "figure3_inputs_verified.txt")
F3_S = script("figure3", "").rstrip("/\\")
F3_PB = "shared/ec_identity_pseudobulk/"
F3_GROUPS = result("figure2", "gene_groups.tsv.gz")


rule figure3_levels:
    input:
        ok=F3_OK,
        groups=F3_GROUPS,
        gut_sums=data(F3_PB + "gut_primary_sums.tsv.gz"),
        gut_groups=data(F3_PB + "gut_primary_groups.tsv"),
        hio_sums=data(F3_PB + "hio_mccauley_endocrine1_sums.tsv.gz"),
        hio_groups=data(F3_PB + "hio_mccauley_endocrine1_groups.tsv"),
        sc_sums=data(F3_PB + "sc_sums.tsv.gz"),
        sc_groups=data(F3_PB + "sc_groups.tsv"),
        islet_sums=data(F3_PB + "islet_sums.tsv.gz"),
        islet_groups=data(F3_PB + "islet_groups.tsv"),
        gut_roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
        islet_roster=data("shared/adult_islet_roster/roster.json"),
    output:
        medians=result("figure3", "group_medians.tsv"),
        g1_levels=result("figure3", "g1_pooled_levels.tsv"),
        g1_counts=result("figure3", "g1_column_counts.tsv"),
    shell:
        '"{RSCRIPT}" "{F3_S}/levels.R" groups="{input.groups}" gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" '
        'hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" sc_sums="{input.sc_sums}" sc_groups="{input.sc_groups}" '
        'islet_sums="{input.islet_sums}" islet_groups="{input.islet_groups}" gut_roster="{input.gut_roster}" '
        'islet_roster="{input.islet_roster}" medians="{output.medians}" g1_levels="{output.g1_levels}" g1_counts="{output.g1_counts}"'


rule figure3_preprint_lists:
    input:
        ok=F3_OK,
        groups=F3_GROUPS,
        liu=data("figure3/preprint_gene_lists/Liu2025_SupplementaryTable5_gene_programs.csv"),
        fiancette=data("figure3/preprint_gene_lists/Fiancette2026_TableS3_cluster_markers.xlsx"),
    output:
        out=result("figure3", "preprint_gene_lists.tsv"),
    shell:
        '"{PYTHON}" "{F3_S}/preprint_lists.py" groups="{input.groups}" liu="{input.liu}" fiancette="{input.fiancette}" '
        'out="{output.out}"'
