# Figure S5: the Figure 2D gene groups with and without their mixed genes; pathway terms (A) and C8 cell-type terms (B).
# Scripts are in workflow/scripts/figure_s5/. The gene groups, mixed genes and tested-gene background are read from the
# deposit (shared/gene_groups/; the same content as the figure2 gene-group results).

FS5_OK = str(RESULTS / "figure_s5_inputs_verified.txt")
FS5_S = script("figure_s5", "").rstrip("/\\")


rule figure_s5_mixed_gene_pathways:
    input:
        ok=FS5_OK,
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        mixed=data("shared/gene_groups/mixed_genes.tsv"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        msigdb=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all=result("figure_s5", "pathways_all_terms.tsv.gz"),
        called=result("figure_s5", "pathways_called.tsv"),
        compare=result("figure_s5", "pathways_group_vs_group_with_mixed.tsv"),
        lists=result("figure_s5", "lists.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{FS5_S}/mixed_gene_pathways.R" groups="{input.groups}" mixed="{input.mixed}" universe="{input.universe}" '
        'msigdb_zip="{input.msigdb}" all="{output.all}" called="{output.called}" compare="{output.compare}" lists="{output.lists}"'


rule figure_s5_mixed_gene_celltypes:
    input:
        ok=FS5_OK,
        lists=result("figure_s5", "lists.tsv"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        msigdb=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all=result("figure_s5", "celltypes_all_terms.tsv"),
        called=result("figure_s5", "celltypes_called.tsv"),
        compare=result("figure_s5", "celltypes_group_vs_group_with_mixed.tsv"),
        by_class=result("figure_s5", "celltypes_called_by_class.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{FS5_S}/mixed_gene_celltypes.R" lists="{input.lists}" universe="{input.universe}" msigdb_zip="{input.msigdb}" '
        'all="{output.all}" called="{output.called}" compare="{output.compare}" by_class="{output.by_class}"'
