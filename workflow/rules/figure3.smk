# Figure 3: tables behind panels A-H. Scripts are in workflow/scripts/figure3/.
# R model fits carry resources heavy_r=1; run with --resources heavy_r=1 to keep them one at a time.

F3_OK = str(RESULTS / "figure3_inputs_verified.txt")
F3_S = script("figure3", "").rstrip("/\\")
F3_AGG = "shared/ec_identity_pseudobulk/"
F3_SC_SAMPLES = [
    "GSE151117__GSE151117__SRR11837751__SC_SIC.tsv.gz",
    "GSE167880__GSE167880__SRR13805895__SC_SIC.tsv.gz",
    "GSE167880__GSE167880__SRR13805896__SC_SIC.tsv.gz",
    "GSE167880__GSE167880__SRR13805897__SC_SIC.tsv.gz",
    "GSE185036_hCiPSC-islets_S6D2__hCiPSC-islets_S6D2__SC_SIC.tsv.gz",
    "GSE202497__GSE202497__SRR19140222__SC_SIC.tsv.gz",
    "GSE202497__GSE202497__SRR19140223__SC_SIC.tsv.gz",
    "GSE226346_MUC18396__S6_WT__SC_SIC.tsv.gz",
    "GSE199636__sc_islet_1.tsv.gz",
    "GSE199636__sc_islet_2.tsv.gz",
    "GSE199636__s6_week2.tsv.gz",
]


def f3_agg(key, stem):
    """The three files of one pseudobulk, as key_sums= key_detected= key_groups= arguments."""
    return {f"{key}_sums": data(F3_AGG + stem + "_sums.tsv.gz"),
            f"{key}_detected": data(F3_AGG + stem + "_detected.tsv.gz"),
            f"{key}_groups": data(F3_AGG + stem + "_groups.tsv")}


def f3_args(named):
    return " ".join(f'{k}="{v}"' for k, v in named.items())


rule figure3_census:
    input:
        ok=F3_OK,
        census=data("figure3/primary_label_population_census.tsv"),
    output:
        out=result("figure3", "panel_a_census.tsv"),
    shell:
        '"{RSCRIPT}" "{F3_S}/census.R" census="{input.census}" out="{output.out}"'


rule figure3_score_by_study:
    input:
        ok=F3_OK,
        scores=data("figure3/f3_03_scores_mature_endocrine.tsv.gz"),
        fit=data("shared/four_marker_classifier/fitted_constants.json"),
        census=data("figure3/primary_label_population_census.tsv"),
    output:
        cells=result("figure3", "panel_b_scores_by_cell.tsv.gz"),
        cutoff=result("figure3", "panel_b_cutoff.tsv"),
        above=result("figure3", "panel_b_above_cutoff.tsv"),
    shell:
        '"{RSCRIPT}" "{F3_S}/score_by_study.R" scores="{input.scores}" fit="{input.fit}" census="{input.census}" '
        'cells="{output.cells}" cutoff="{output.cutoff}" above="{output.above}"'


F3_C_IN = {**f3_agg("gut", "gut_primary"), **f3_agg("sc", "sc"), **f3_agg("ind", "induced"),
           **f3_agg("hio", "hio_mccauley_endocrine1"), "roster": data("shared/gut_reference_roster/gut_reference_roster.json")}

rule figure3_signature_correlations:
    input:
        ok=F3_OK,
        **F3_C_IN,
    output:
        dots=result("figure3", "panel_c_dots.tsv"),
        hio_dots=result("figure3", "panel_c_hio_dots.tsv"),
        pairwise=result("figure3", "panel_c_pairwise.tsv"),
        hio_pairwise=result("figure3", "panel_c_hio_pairwise.tsv"),
        gut_folds=result("figure3", "gut_ec_fold_changes.tsv.gz"),
    params:
        args=f3_args(F3_C_IN),
    shell:
        '"{PYTHON}" "{F3_S}/signature_correlations.py" {params.args} dots="{output.dots}" hio_dots="{output.hio_dots}" '
        'pairwise="{output.pairwise}" hio_pairwise="{output.hio_pairwise}" gut_folds="{output.gut_folds}"'


rule figure3_donor_fit:
    input:
        ok=F3_OK,
        gut_sums=data(F3_AGG + "gut_primary_sums.tsv.gz"),
        gut_groups=data(F3_AGG + "gut_primary_groups.tsv"),
        roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
        samples=[data("shared/sc_sic_sample_pseudobulk/" + f) for f in F3_SC_SAMPLES],
    output:
        levels=result("figure3", "panel_d_donor_fit.tsv.gz"),
        per_study=result("figure3", "panel_e_per_study_levels.tsv.gz"),
        profiles=result("figure3", "pooled_voom_profiles.tsv.gz"),
        tested=result("figure3", "tested_genes.txt"),
    params:
        sc_dir=data("shared/sc_sic_sample_pseudobulk"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{F3_S}/donor_fit.R" gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" '
        'roster="{input.roster}" sc_dir="{params.sc_dir}" levels="{output.levels}" per_study="{output.per_study}" '
        'profiles="{output.profiles}" tested="{output.tested}"'


F3_F_IN = {**f3_agg("sc", "sc"), "gut_sums": data(F3_AGG + "gut_primary_sums.tsv.gz"),
           "gut_groups": data(F3_AGG + "gut_primary_groups.tsv"), "roster": data("shared/gut_reference_roster/gut_reference_roster.json")}

rule figure3_gene_groups:
    input:
        ok=F3_OK,
        gut_folds=result("figure3", "gut_ec_fold_changes.tsv.gz"),
        fit=result("figure3", "panel_d_donor_fit.tsv.gz"),
        profiles=result("figure3", "pooled_voom_profiles.tsv.gz"),
        **F3_F_IN,
    output:
        definitions=result("figure3", "panel_f_gene_groups.tsv.gz"),
        mixed=result("figure3", "panel_f_mixed_genes.tsv"),
        universe=result("figure3", "panel_f_universe.txt"),
    params:
        args=f3_args(F3_F_IN),
    shell:
        '"{PYTHON}" "{F3_S}/gene_groups.py" {params.args} gut_folds="{input.gut_folds}" fit="{input.fit}" '
        'profiles="{input.profiles}" definitions="{output.definitions}" mixed="{output.mixed}" universe="{output.universe}"'


rule figure3_enrich_pathways:
    input:
        ok=F3_OK,
        groups=result("figure3", "panel_f_gene_groups.tsv.gz"),
        universe=result("figure3", "panel_f_universe.txt"),
        fit=result("figure3", "panel_d_donor_fit.tsv.gz"),
        msigdb=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all_terms=result("figure3", "panel_g_pathways_all_terms.tsv.gz"),
        called=result("figure3", "panel_g_pathways_called.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{F3_S}/enrich_pathways.R" groups="{input.groups}" universe="{input.universe}" '
        'fit="{input.fit}" msigdb_zip="{input.msigdb}" all_terms="{output.all_terms}" called="{output.called}"'


rule figure3_enrich_c8:
    input:
        ok=F3_OK,
        groups=result("figure3", "panel_f_gene_groups.tsv.gz"),
        universe=result("figure3", "panel_f_universe.txt"),
        fit=result("figure3", "panel_d_donor_fit.tsv.gz"),
        msigdb=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all_terms=result("figure3", "panel_g_celltypes_all_terms.tsv"),
        called=result("figure3", "panel_g_celltypes_called.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{F3_S}/enrich_c8.R" groups="{input.groups}" universe="{input.universe}" '
        'fit="{input.fit}" msigdb_zip="{input.msigdb}" all_terms="{output.all_terms}" called="{output.called}"'


rule figure3_age_genes:
    input:
        ok=F3_OK,
        gut_sums=data(F3_AGG + "gut_primary_sums.tsv.gz"),
        gut_detected=data(F3_AGG + "gut_primary_detected.tsv.gz"),
        gut_groups=data(F3_AGG + "gut_primary_groups.tsv"),
        excluded=data("figure3/q1_excluded_marker_genes.tsv"),
    output:
        adult=result("figure3", "age_adult_high_genes.txt"),
        fetal=result("figure3", "age_fetal_high_genes.txt"),
        adult_no_floor=result("figure3", "age_adult_high_genes_no_floor.txt"),
        fetal_no_floor=result("figure3", "age_fetal_high_genes_no_floor.txt"),
        calls=result("figure3", "age_calls.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{F3_S}/age_genes.R" gut_sums="{input.gut_sums}" gut_detected="{input.gut_detected}" '
        'gut_groups="{input.gut_groups}" excluded="{input.excluded}" adult="{output.adult}" fetal="{output.fetal}" '
        'adult_no_floor="{output.adult_no_floor}" fetal_no_floor="{output.fetal_no_floor}" calls="{output.calls}"'


F3_H_TEST = "figure3/h_test_sources_four_marker/four_marker_abstain"
F3_H_IN = {f"{key}_{part}": data(F3_AGG + stem + ("_groups.tsv" if part == "groups" else "_sums.tsv.gz"))
           for key, stem in (("gut", "gut_primary"), ("sc", "sc"), ("ind", "induced"), ("hio", "hio_mccauley_endocrine1"))
           for part in ("sums", "groups")}
F3_H_IN.update({"test_sums": data(F3_H_TEST + "_sums.tsv.gz"), "test_groups": data(F3_H_TEST + "_groups.tsv")})


rule figure3_age_scores:
    input:
        ok=F3_OK,
        adult=result("figure3", "age_adult_high_genes.txt"),
        fetal=result("figure3", "age_fetal_high_genes.txt"),
        **F3_H_IN,
    output:
        drawn=result("figure3", "panel_h_drawn.tsv"),
        genes=result("figure3", "panel_h_genes.tsv"),
    params:
        args=f3_args(F3_H_IN),
    shell:
        '"{PYTHON}" "{F3_S}/age_scores.py" {params.args} adult="{input.adult}" fetal="{input.fetal}" '
        'drawn="{output.drawn}" genes="{output.genes}"'
