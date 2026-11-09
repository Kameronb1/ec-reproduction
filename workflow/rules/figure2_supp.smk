# Figure 2 supplement: the four-marker EC classifier's fit in the fitting studies (S1), its maps
# and clusters (S2), the cells it classifies as EC in the test studies and HIO grafts (S3), the 20-gene
# classifier's gene selection and controls (S4), and how it classifies a cell, its weights and both
# classifiers per donor or sample (S5). Figure 2's leave-one-study-out run of the 20-gene
# classifier is a deposited intermediate (figure2/classifier_run/ and figure2_supp/classifier_run/).
# The S4C matched-random test (200 random gene sets refitted per study, about 55 min on one core) is
# deposited too and is made by classifier.R when given matched_random= and hio_matched= outputs;
# everything else runs locally.
# Rules marked resources heavy_r=1 hold the fitting-cell matrix, a test study's full count matrix
# or the training corpus in memory; with --resources heavy_r=1 they run one at a time. Their two
# threads are passed by Snakemake to the BLAS and OpenMP libraries the R packages use.

F2S_OK = str(RESULTS / "figure2_supp_inputs_verified.txt")
F2S_MATRIX = dict(
    matrix=data("shared/fitting_matrix/fitting_matrix.f64le"),
    features=data("shared/fitting_matrix/features.tsv"),
    cells=data("shared/fitting_matrix/cells.tsv"),
    meta=data("shared/fitting_matrix/cell_metadata.tsv.gz"),
)
F2S_H5AD = [data(f"figure2/heldout_studies/{s}.h5ad") for s in ("gut_cell_atlas", "hickey", "egozi", "gse238274")]
F2S_H5AD_DIR = str(Path(F2S_H5AD[0]).parent)
F2S_ANALYSIS = data("shared/gut_reference_analysis.rds")
F2S_HIO = dict(
    hio_counts=data("shared/gse214852_hio/counts.mtx"),
    hio_genes=data("shared/gse214852_hio/genes.txt"),
    hio_cells=data("shared/gse214852_hio/cells.tsv"),
)


def f2s(name):
    return result("figure2_supp", name)


rule figure2_supp_marker_fit:
    input:
        ok=F2S_OK,
        **F2S_MATRIX,
    output:
        cells_out=f2s("fitting_cell_classes.tsv.gz"),
        counts_out=f2s("fitting_class_counts.tsv"),
        arms_out=f2s("fitting_arm_summary.tsv"),
        excl_out=f2s("fitting_exclusion_set_comparison.tsv"),
        fit_out=f2s("intermediate/four_marker_fit.json"),
        hio_cells_out=f2s("hpsc_cell_classes.tsv.gz"),
        hio_arms_out=f2s("hpsc_arm_summary.tsv"),
        summary_out=f2s("score_summary_by_arm_with_hpsc.tsv"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "marker_fit.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" matrix="{input.matrix}" features="{input.features}" cells="{input.cells}" '
        'meta="{input.meta}" cells_out="{output.cells_out}" counts_out="{output.counts_out}" '
        'arms_out="{output.arms_out}" excl_out="{output.excl_out}" fit_out="{output.fit_out}" '
        'hio_cells_out="{output.hio_cells_out}" hio_arms_out="{output.hio_arms_out}" summary_out="{output.summary_out}"'


rule figure2_supp_outcome_by_group:
    input:
        ok=F2S_OK,
        cells=f2s("fitting_cell_classes.tsv.gz"),
    output:
        f2s("fitting_outcome_by_arm.tsv"),
    params:
        script=script("figure2_supp", "outcome_by_group.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" out="{output}"'


rule figure2_supp_cutoff_stability:
    input:
        ok=F2S_OK,
        **F2S_MATRIX,
    output:
        f2s("cutoff_stability_summary.tsv"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "cutoff_stability.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" matrix="{input.matrix}" features="{input.features}" cells="{input.cells}" '
        'meta="{input.meta}" out="{output}"'


rule figure2_supp_umap:
    input:
        ok=F2S_OK,
        classes=f2s("fitting_cell_classes.tsv.gz"),
        **F2S_MATRIX,
    output:
        coords=f2s("umap_coordinates.tsv.gz"),
        purity=f2s("umap_knn_purity.tsv"),
        scaling=f2s("panel_gene_scaling.tsv"),
        embeddings=f2s("harmony_embeddings.tsv.gz"),
        gene_expression=f2s("gene_expression.tsv"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "umap.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" matrix="{input.matrix}" features="{input.features}" cells="{input.cells}" '
        'meta="{input.meta}" classes="{input.classes}" coords="{output.coords}" purity="{output.purity}" '
        'scaling="{output.scaling}" embeddings="{output.embeddings}" gene_expression="{output.gene_expression}"'


rule figure2_supp_clusters:
    input:
        ok=F2S_OK,
        classes=f2s("fitting_cell_classes.tsv.gz"),
        embeddings=f2s("harmony_embeddings.tsv.gz"),
        **F2S_MATRIX,
    output:
        sweep=f2s("resolution_sweep.tsv"),
        summary=f2s("cluster_summary.tsv"),
        markers=f2s("cluster_markers.tsv.gz"),
        assign=f2s("cell_clusters.tsv.gz"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "clusters.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" matrix="{input.matrix}" features="{input.features}" cells="{input.cells}" '
        'meta="{input.meta}" classes="{input.classes}" embeddings="{input.embeddings}" sweep="{output.sweep}" '
        'summary="{output.summary}" markers="{output.markers}" assign="{output.assign}"'


rule figure2_supp_heldout:
    input:
        ok=F2S_OK,
        h5ad=F2S_H5AD,
        fit=f2s("intermediate/four_marker_fit.json"),
        stability=f2s("cutoff_stability_summary.tsv"),
        meta=data("shared/fitting_matrix/cell_metadata.tsv.gz"),
    output:
        scoring=f2s("intermediate/heldout_scoring.tsv.gz"),
        by_arm=f2s("intermediate/scoring_by_arm.tsv"),
        coords=f2s("heldout_space_coords.tsv.gz"),
        settings=f2s("heldout_space_settings.tsv"),
        sweep=f2s("cutoff_sweep.tsv.gz"),
        landmarks=f2s("sweep_landmarks.tsv"),
        refits=f2s("transported_refits.tsv"),
        overlap=f2s("refit_cell_overlap.tsv"),
        comparisons=f2s("refit_cell_comparisons.tsv.gz"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "heldout.R"),
        h5ad_dir=F2S_H5AD_DIR,
    shell:
        '"{RSCRIPT}" "{params.script}" h5ad_dir="{params.h5ad_dir}" fit="{input.fit}" stability="{input.stability}" '
        'meta="{input.meta}" scoring="{output.scoring}" by_arm="{output.by_arm}" coords="{output.coords}" '
        'settings="{output.settings}" sweep="{output.sweep}" landmarks="{output.landmarks}" '
        'refits="{output.refits}" overlap="{output.overlap}" comparisons="{output.comparisons}"'


rule figure2_supp_hio_four_marker:
    input:
        ok=F2S_OK,
        fit=f2s("intermediate/four_marker_fit.json"),
        stability=f2s("cutoff_stability_summary.tsv"),
        hio_classes=f2s("hpsc_cell_classes.tsv.gz"),
        hio_cells=data("shared/gse214852_hio/cells.tsv"),
        refits=f2s("transported_refits.tsv"),
        landmarks=f2s("sweep_landmarks.tsv"),
        **F2S_MATRIX,
    output:
        coords=f2s("hio_heldout_space_coords.tsv.gz"),
        sweep=f2s("hio_cutoff_sweep.tsv.gz"),
        refits_out=f2s("hio_transported_refits.tsv"),
        landmarks_out=f2s("hio_sweep_landmarks.tsv"),
        depth=f2s("hio_four_marker_depth.tsv"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "hio_four_marker.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" matrix="{input.matrix}" features="{input.features}" cells="{input.cells}" '
        'meta="{input.meta}" fit="{input.fit}" stability="{input.stability}" hio_classes="{input.hio_classes}" '
        'hio_cells="{input.hio_cells}" refits="{input.refits}" landmarks="{input.landmarks}" '
        'coords="{output.coords}" sweep="{output.sweep}" refits_out="{output.refits_out}" '
        'landmarks_out="{output.landmarks_out}" depth="{output.depth}"'


rule figure2_supp_four_marker_depth:
    input:
        ok=F2S_OK,
        h5ad=F2S_H5AD,
        fit=f2s("intermediate/four_marker_fit.json"),
        scoring=f2s("intermediate/heldout_scoring.tsv.gz"),
    output:
        f2s("four_marker_depth.tsv"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "four_marker_depth.R"),
        h5ad_dir=F2S_H5AD_DIR,
    shell:
        '"{RSCRIPT}" "{params.script}" h5ad_dir="{params.h5ad_dir}" fit="{input.fit}" scoring="{input.scoring}" '
        'out="{output}"'


rule figure2_supp_classifier:
    input:
        ok=F2S_OK,
        core=data("figure2/classifier_run/candidate_panel_core_result.rds"),
        analysis=F2S_ANALYSIS,
        secondary=data("figure2/classifier_inputs/secondary_validation.rds"),
        singh_labels=data("figure2/classifier_run/sv_labels.tsv.gz"),
        hio_classes=f2s("hpsc_cell_classes.tsv.gz"),
        **F2S_HIO,
    output:
        predictions=f2s("predictions.tsv.gz"),
        weights=f2s("weights.tsv"),
        hio_predictions=f2s("hio_classifier_predictions.tsv.gz"),
        hio_weights=f2s("hio_classifier_weights.tsv"),
        generic=f2s("generic_neuroendocrine.tsv"),
        hio_generic=f2s("hio_generic_neuroendocrine.tsv"),
        depth=f2s("depth.tsv"),
        hio_depth=f2s("hio_classifier_depth.tsv"),
        hio_by_organoid=f2s("hio_calls_by_organoid.tsv"),
        hio_ranking=f2s("hio_ranking_versus_calling.tsv"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "classifier.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" core="{input.core}" analysis="{input.analysis}" secondary="{input.secondary}" '
        'singh_labels="{input.singh_labels}" hio_counts="{input.hio_counts}" hio_genes="{input.hio_genes}" '
        'hio_cells="{input.hio_cells}" hio_classes="{input.hio_classes}" predictions="{output.predictions}" '
        'weights="{output.weights}" hio_predictions="{output.hio_predictions}" hio_weights="{output.hio_weights}" '
        'generic="{output.generic}" '
        'hio_generic="{output.hio_generic}" depth="{output.depth}" hio_depth="{output.hio_depth}" '
        'hio_by_organoid="{output.hio_by_organoid}" hio_ranking="{output.hio_ranking}"'


rule figure2_supp_matched_random:
    """The S4C matched-random tables, deposited (made by classifier.R with matched_random= and hio_matched=),
    copied into the results as they are so the plot reads them like every other table."""
    input:
        ok=F2S_OK,
        gut=data("figure2_supp/classifier_run/matched_random.tsv"),
        hio=data("figure2_supp/classifier_run/hio_matched_random.tsv"),
    output:
        gut=f2s("matched_random.tsv"),
        hio=f2s("hio_matched_random.tsv"),
    shell:
        '"{PYTHON}" -c "import shutil, sys; shutil.copyfile(sys.argv[1], sys.argv[2]); shutil.copyfile(sys.argv[3], sys.argv[4])" '
        '"{input.gut}" "{output.gut}" "{input.hio}" "{output.hio}"'


rule figure2_supp_gene_selection:
    input:
        ok=F2S_OK,
        panels=data("figure2/classifier_run/nested_panels.tsv.gz"),
        d_genes=data("shared/figure2_panel_d/D_effect_heatmap.tsv"),
        replicates=data("figure2_supp/classifier_run/recurrence_null_replicates.tsv"),
    output:
        grid=f2s("gene_selection_by_fold.tsv"),
        null=f2s("recurrence_null.tsv"),
        observed=f2s("observed_recurrence.tsv"),
        summary=f2s("recurrence_summary.tsv"),
    params:
        script=script("figure2_supp", "gene_selection.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" panels="{input.panels}" d_genes="{input.d_genes}" '
        'replicates="{input.replicates}" grid="{output.grid}" null="{output.null}" observed="{output.observed}" '
        'summary="{output.summary}"'


rule figure2_supp_classifier_explainer:
    input:
        ok=F2S_OK,
        analysis=F2S_ANALYSIS,
        weights=f2s("weights.tsv"),
        predictions=f2s("predictions.tsv.gz"),
        d_genes=data("shared/figure2_panel_d/D_effect_heatmap.tsv"),
    output:
        genes=f2s("classifier_explainer_genes.tsv"),
        mixture=f2s("classifier_explainer_mixture.tsv"),
        hist=f2s("classifier_explainer_score_histogram.tsv"),
        mix=f2s("classifier_explainer_mixture_parameters.tsv"),
    threads: 2
    resources:
        heavy_r=1,
    params:
        script=script("figure2_supp", "classifier_explainer.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" analysis="{input.analysis}" weights="{input.weights}" '
        'predictions="{input.predictions}" d_genes="{input.d_genes}" genes="{output.genes}" mixture="{output.mixture}" '
        'hist="{output.hist}" mix="{output.mix}"'


rule figure2_supp_per_sample:
    input:
        ok=F2S_OK,
        scoring=f2s("intermediate/heldout_scoring.tsv.gz"),
        predictions=f2s("predictions.tsv.gz"),
    output:
        f2s("calls_by_donor.tsv"),
    params:
        script=script("figure2_supp", "per_sample.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" scoring="{input.scoring}" predictions="{input.predictions}" out="{output}"'


# The TPH1-only classifier (S5C, S6A): Figure 2's four-marker fit and test-study scoring rerun with TPH1 alone
# (markers=TPH1 min_detected=1), into tph1_only/. Its histogram and mixture tables are not drawn.
rule figure2_supp_tph1_fit:
    input:
        ok=F2S_OK,
        **F2S_MATRIX,
    output:
        hist=f2s("tph1_only/A_histogram.tsv"),
        curves=f2s("tph1_only/A_mixture_curves.tsv"),
        fit=f2s("tph1_only/four_marker_fit.json"),
        hio=f2s("tph1_only/hio_four_marker_calls.tsv"),
    params:
        script=script("figure2", "four_marker_fit.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" matrix="{input.matrix}" features="{input.features}" '
        'cells="{input.cells}" meta="{input.meta}" hist="{output.hist}" curves="{output.curves}" '
        'fit="{output.fit}" hio="{output.hio}" markers=TPH1 min_detected=1'


rule figure2_supp_tph1_heldout:
    input:
        ok=F2S_OK,
        fit=f2s("tph1_only/four_marker_fit.json"),
        h5ad=F2S_H5AD,
    output:
        units=f2s("tph1_only/B_units.tsv"),
        labels=f2s("tph1_only/C_primary_labels.tsv"),
        cells=f2s("tph1_only/heldout_cells.tsv"),
    params:
        script=script("figure2", "heldout_marker_calls.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" fit="{input.fit}" '
        'h5ad_gut_cell_atlas="{input.h5ad[0]}" h5ad_hickey="{input.h5ad[1]}" h5ad_egozi="{input.h5ad[2]}" '
        'h5ad_gse238274="{input.h5ad[3]}" units="{output.units}" labels="{output.labels}" cells="{output.cells}"'


rule figure2_supp_tph1_calls:
    input:
        ok=F2S_OK,
        cells4=result("figure2", "intermediate/heldout_cells.tsv"),
        cells1=f2s("tph1_only/heldout_cells.tsv"),
        hio4=result("figure2", "intermediate/hio_four_marker_calls.tsv"),
        hio1=f2s("tph1_only/hio_four_marker_calls.tsv"),
        predictions=f2s("predictions.tsv.gz"),
        hio_predictions=f2s("hio_classifier_predictions.tsv.gz"),
        by_donor=f2s("calls_by_donor.tsv"),
        by_organoid=f2s("hio_calls_by_organoid.tsv"),
    output:
        units=f2s("tph1_calls_by_unit.tsv"),
        other=f2s("other_eec_calls.tsv"),
    params:
        script=script("figure2_supp", "tph1_calls.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells4="{input.cells4}" cells1="{input.cells1}" hio4="{input.hio4}" '
        'hio1="{input.hio1}" predictions="{input.predictions}" hio_predictions="{input.hio_predictions}" '
        'by_donor="{input.by_donor}" by_organoid="{input.by_organoid}" units="{output.units}" other="{output.other}"'
