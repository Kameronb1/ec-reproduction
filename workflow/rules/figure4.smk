# Figure 4 and its supplementary Figures S14 and S16: SC-islet composition (B, S14), per-cell EC scores and expression in grafts
# and longer culture (C), cholesterol synthesis (D, S16A), gap to primary intestinal EC (E), fetal-versus-adult age score (F),
# NFE2L1 regulon target score (S16A) and regulon screens (S16B-D). Panel A is the curated design of the rows (samples in
# rows_called.tsv).

FIG4_OK = str(RESULTS / "figure4_inputs_verified.txt")
F4_TREES = ["GSE167880", "GSE151117", "GSE139535", "GSE185036_S6D2", "GSE185036_10wpt", "GSE202497", "GSE114412"]


def f4_tree_files(n):
    """Inputs tree1..tree<n> (profiles, genes, counts) in the order the scripts expect."""
    files = {}
    for i, t in enumerate(F4_TREES[:n], start=1):
        for part, name in [("profiles", "profiles.tsv"), ("genes", "genes.tsv"), ("counts", "counts.tsv.gz")]:
            files[f"tree{i}_{part}"] = data(f"figure4/pseudobulk/{t}/{name}")
    return files


def f4_args(wildcards, input):
    """name="path" for every named input except the input-check file (ok)."""
    return " ".join(f'{k}="{v}"' for k, v in input.items() if k != "ok")


rule figure4_composition:
    input:
        ok=FIG4_OK,
        balboa=data("figure4/pseudobulk/GSE167880/profiles.tsv"),
        hues8=data("figure4/pseudobulk/GSE151117/profiles.tsv"),
        ws4=data("figure4/pseudobulk/GSE139535/profiles.tsv"),
        du_s6d2=data("figure4/pseudobulk/GSE185036_S6D2/profiles.tsv"),
        du_10wpt=data("figure4/pseudobulk/GSE185036_10wpt/profiles.tsv"),
        zhu=data("figure4/zhu2023_primary_labels.tsv.gz"),
        a23_profiles=data("figure4/gse199636_recount/profiles.tsv"),
        a23_groups=data("figure4/gse199636_primary_label_groups.tsv"),
        extra=data("figure4/abundance_extra_studies.tsv"),
    output:
        f4b=result("figure4", "f4b_composition.tsv"),
        s1=result("figure4", "s1_composition.tsv"),
    params:
        script=script("figure4", "composition.R"),
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{params.script}" {params.args} f4b="{output.f4b}" s1="{output.s1}"'


rule figure4_ec_calls:
    input:
        ok=FIG4_OK,
        calls=data("figure4/ec_calls_per_cell.tsv.gz"),
        calls_soupx=data("figure4/ec_calls_per_cell_soupx.tsv.gz"),
        comparisons=data("figure4/comparisons_whole_cell.tsv"),
        comparisons_nuclei=data("figure4/comparisons_nuclei.tsv"),
    output:
        rows=result("figure4", "rows_called.tsv"),
        rows_soupx=result("figure4", "rows_called_soupx.tsv"),
        cells=result("figure4", "C_cells.tsv.gz"),
    params:
        script=script("figure4", "ec_calls.R"),
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{params.script}" {params.args} rows="{output.rows}" rows_soupx="{output.rows_soupx}" cells="{output.cells}"'


rule figure4_levels:
    input:
        ok=FIG4_OK,
        **f4_tree_files(6),
        gut=data("shared/gut_reference_analysis.rds"),
        cfap92=data("figure4/gut_cfap92_alias_row.rds"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        nuc_profiles=data("figure4/gse199636_recount/profiles.tsv"),
        nuc_pseudobulk=data("figure4/gse199636_recount/pseudobulk_counts.tsv.gz"),
        soupx_profiles=data("figure4/soupx_graft_profiles.tsv"),
        soupx_pseudobulk=data("figure4/soupx_graft_pseudobulk.tsv.gz"),
        rows=result("figure4", "rows_called.tsv"),
        classifier=data("figure4/sc_sic_classifier_fold_models.rds"),
        chol=data("figure4/cholesterol_biosynthesis_genes.tsv"),
        groups=data("shared/fig3_gene_groups.tsv.gz"),
        roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
    output:
        bars=result("figure4", "C_bars.tsv"),
        chol_out=result("figure4", "D_cholesterol_medians.tsv"),
        chol_genes=result("figure4", "D_cholesterol_per_gene.tsv"),
        gaps=result("figure4", "E_group_gaps.tsv"),
        end_spread=result("figure4", "E_end_spread.tsv"),
    params:
        script=script("figure4", "levels.R"),
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{params.script}" {params.args} bars="{output.bars}" chol_out="{output.chol_out}" '
        'chol_genes="{output.chol_genes}" gaps="{output.gaps}" end_spread="{output.end_spread}"'


F4_AGE_LISTS = {"adult": data("figure4/age_adult_high_genes.txt"), "fetal": data("figure4/age_fetal_high_genes.txt")}
F4_AGE_FLOOR = {"adult_floor": result("figure4", "age_adult_high_genes_floor10.txt"),
                "fetal_floor": result("figure4", "age_fetal_high_genes_floor10.txt")}
F4_GUT_PB = "shared/ec_identity_pseudobulk/gut_primary"
F4_TEST_PB = "figure3/h_test_sources_four_marker/four_marker_abstain"


rule figure4_age_floor_genes:
    input:
        ok=FIG4_OK,
        **F4_AGE_LISTS,
        gut_detected=data(F4_GUT_PB + "_detected.tsv.gz"),
        gut_groups=data(F4_GUT_PB + "_groups.tsv"),
    output:
        **F4_AGE_FLOOR,
    params:
        script=script("figure4", "age_floor_genes.R"),
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{params.script}" {params.args} adult_floor="{output.adult_floor}" fetal_floor="{output.fetal_floor}"'


rule figure4_age_score:
    input:
        ok=FIG4_OK,
        **f4_tree_files(5),
        nuc_profiles=data("figure4/gse199636_recount/profiles.tsv"),
        nuc_pseudobulk=data("figure4/gse199636_recount/pseudobulk_counts.tsv.gz"),
        gut_sums=data(F4_GUT_PB + "_sums.tsv.gz"),
        gut_groups=data(F4_GUT_PB + "_groups.tsv"),
        test_sums=data(F4_TEST_PB + "_sums.tsv.gz"),
        test_groups=data(F4_TEST_PB + "_groups.tsv"),
        **F4_AGE_FLOOR,
        rows=result("figure4", "rows_called.tsv"),
    output:
        out=result("figure4", "age_score_by_sample.tsv"),
        ticks=result("figure4", "fig3H_ticks.tsv"),
        genes=result("figure4", "age_score_genes_by_sample.tsv"),
    params:
        script=script("figure4", "age_score.R"),
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{params.script}" {params.args} out="{output.out}" ticks="{output.ticks}" genes="{output.genes}"'


rule figure4_cholesterol_culture:
    input:
        ok=FIG4_OK,
        **f4_tree_files(7),
        gut=data("shared/gut_reference_analysis.rds"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        nuc_profiles=data("figure4/gse199636_recount/profiles.tsv"),
        nuc_pseudobulk=data("figure4/gse199636_recount/pseudobulk_counts.tsv.gz"),
        soupx_profiles=data("figure4/soupx_graft_profiles.tsv"),
        soupx_pseudobulk=data("figure4/soupx_graft_pseudobulk.tsv.gz"),
        rows=result("figure4", "rows_called.tsv"),
        comparisons=data("figure4/comparisons_whole_cell.tsv"),
        chol=data("figure4/cholesterol_biosynthesis_genes.tsv"),
        zhu_alpha_profiles=data("figure4/sc_alpha_profiles/GSE202497/profiles.tsv"),
        zhu_alpha_counts=data("figure4/sc_alpha_profiles/GSE202497/counts.tsv.gz"),
        veres_alpha_profiles=data("figure4/sc_alpha_profiles/GSE114412/profiles.tsv"),
        veres_alpha_counts=data("figure4/sc_alpha_profiles/GSE114412/counts.tsv.gz"),
        d_genes=result("figure4", "D_cholesterol_per_gene.tsv"),
        d_medians=result("figure4", "D_cholesterol_medians.tsv"),
    output:
        tiles=result("figure4", "s2_cholesterol_tiles.tsv"),
        medians=result("figure4", "s2_cholesterol_medians.tsv"),
    params:
        script=script("figure4", "cholesterol_culture.R"),
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{params.script}" {params.args} tiles="{output.tiles}" medians="{output.medians}"'


rule figure4_regulon_screen_tables:
    input:
        ok=FIG4_OK,
        screen=data("figure4/regulon_screen.tsv.gz"),
    output:
        a=result("figure4", "s3a_decline_screen.tsv"),
        b=result("figure4", "s3b_cholesterol_screen.tsv"),
        d=result("figure4", "s3d_regulon_recurrence.tsv"),
    params:
        script=script("figure4", "regulon_screen_tables.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" screen="{input.screen}" a="{output.a}" b="{output.b}" d="{output.d}"'


rule figure4_nfe2l1_activity:
    input:
        ok=FIG4_OK,
        units=data("figure4/nfe2l1_regulon_unit_activity.tsv"),
    output:
        result("figure4", "s2_nfe2l1_bars.tsv"),
    params:
        script=script("figure4", "nfe2l1_activity.py"),
    shell:
        '"{PYTHON}" "{params.script}" units="{input.units}" out="{output}"'
