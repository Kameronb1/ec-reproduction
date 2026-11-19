# Supplementary Figure S8: cholesterol synthesis genes in grafts and extended culture (A). Panels B and C draw the
# chromatin tables of the Figure 4 unit; this rule reads that unit's row table (samples at each row end).
# The rule holds the intestinal reference matrix in memory (resources heavy_r=1).

FIGS8_OK = str(RESULTS / "figure_s8_inputs_verified.txt")
S8_TREES = ["GSE167880", "GSE151117", "GSE139535", "GSE185036_S6D2", "GSE185036_10wpt", "GSE202497", "GSE114412"]


def s8_tree_files():
    """Inputs tree1..tree7 (profiles, genes, counts) in the order the script expects."""
    files = {}
    for i, t in enumerate(S8_TREES, start=1):
        for part, name in [("profiles", "profiles.tsv"), ("genes", "genes.tsv"), ("counts", "counts.tsv.gz")]:
            files[f"tree{i}_{part}"] = data(f"figure4/pseudobulk/{t}/{name}")
    return files


rule figure_s8_cholesterol_culture:
    """Figure S8A: change of the 26 measured Reactome cholesterol biosynthesis genes per row and cell type, and their
    median."""
    input:
        ok=FIGS8_OK,
        script=script("figure_s8", "cholesterol_culture.R"),
        **s8_tree_files(),
        gut=data("shared/gut_reference_analysis.rds"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        nuc_profiles=data("figure4/gse199636_recount/profiles.tsv"),
        nuc_pseudobulk=data("figure4/gse199636_recount/pseudobulk_counts.tsv.gz"),
        soupx_profiles=data("figure4/soupx_graft_profiles.tsv"),
        soupx_pseudobulk=data("figure4/soupx_graft_pseudobulk.tsv.gz"),
        rows=result("figure4", "row_samples.tsv"),
        comparisons=data("figure4/comparisons_whole_cell.tsv"),
        chol=data("figure4/cholesterol_biosynthesis_genes.tsv"),
        zhu_alpha_profiles=data("figure4/sc_alpha_profiles/GSE202497/profiles.tsv"),
        zhu_alpha_counts=data("figure4/sc_alpha_profiles/GSE202497/counts.tsv.gz"),
        veres_alpha_profiles=data("figure4/sc_alpha_profiles/GSE114412/profiles.tsv"),
        veres_alpha_counts=data("figure4/sc_alpha_profiles/GSE114412/counts.tsv.gz"),
    output:
        tiles=result("figure_s8", "s2_cholesterol_tiles.tsv"),
        medians=result("figure_s8", "s2_cholesterol_medians.tsv"),
    params:
        args=lambda wildcards, input: " ".join(f'{k}="{v}"' for k, v in input.items() if k not in ("ok", "script")),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{input.script}" {params.args} tiles="{output.tiles}" medians="{output.medians}"'
