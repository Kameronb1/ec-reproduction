# Figure 2: fit of the label-free four-marker EC classifier (A), the cells it classifies as EC across
# donors and samples (B) and within primary study labels, with the four markers' detection per label (C),
# the gene-ranking datasets' 20-gene panel, its per-study effects and detection (D), and the label-trained
# 20-gene classifier against the four-marker classifier, cell by cell (E).
# The leave-one-study-out run of the 20-gene classifier was run on a computing cluster; its results are
# deposited (figure2/classifier_run/ under data_dir) and read as inputs. Everything else runs locally.

FIG2_OK = str(RESULTS / "figure2_inputs_verified.txt")
FIG2_RUN = "figure2/classifier_run/"


rule figure2_four_marker_fit:
    input:
        ok=FIG2_OK,
        matrix=data("shared/fitting_matrix/fitting_matrix.f64le"),
        features=data("shared/fitting_matrix/features.tsv"),
        cells=data("shared/fitting_matrix/cells.tsv"),
        meta=data("shared/fitting_matrix/cell_metadata.tsv.gz"),
    output:
        hist=result("figure2", "A_histogram.tsv"),
        curves=result("figure2", "A_mixture_curves.tsv"),
        fit=result("figure2", "intermediate/four_marker_fit.json"),
        hio=result("figure2", "intermediate/hio_four_marker_calls.tsv"),
    params:
        script=script("figure2", "four_marker_fit.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" matrix="{input.matrix}" features="{input.features}" '
        'cells="{input.cells}" meta="{input.meta}" hist="{output.hist}" curves="{output.curves}" '
        'fit="{output.fit}" hio="{output.hio}"'


rule figure2_heldout_marker_calls:
    input:
        ok=FIG2_OK,
        fit=result("figure2", "intermediate/four_marker_fit.json"),
        gca=data("figure2/heldout_studies/gut_cell_atlas.h5ad"),
        hickey=data("figure2/heldout_studies/hickey.h5ad"),
        egozi=data("figure2/heldout_studies/egozi.h5ad"),
        singh=data("figure2/heldout_studies/gse238274.h5ad"),
    output:
        units=result("figure2", "B_units.tsv"),
        labels=result("figure2", "C_primary_labels.tsv"),
        cells=result("figure2", "intermediate/heldout_cells.tsv"),
    params:
        script=script("figure2", "heldout_marker_calls.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" fit="{input.fit}" '
        'h5ad_gut_cell_atlas="{input.gca}" h5ad_hickey="{input.hickey}" h5ad_egozi="{input.egozi}" '
        'h5ad_gse238274="{input.singh}" units="{output.units}" labels="{output.labels}" cells="{output.cells}"'


rule figure2_gene_panel_effects:
    input:
        ok=FIG2_OK,
        analysis=data("shared/gut_reference_analysis.rds"),
        secondary=data("figure2/classifier_inputs/secondary_validation.rds"),
        run_metadata=data(FIG2_RUN + "training_metadata.tsv.gz"),
        nested_predictions=data(FIG2_RUN + "nested_predictions.tsv.gz"),
        sv_predictions=data(FIG2_RUN + "sv_predictions.tsv.gz"),
        sv_labels=data(FIG2_RUN + "sv_labels.tsv.gz"),
        sv_panels=data(FIG2_RUN + "sv_panels.tsv"),
        labels=result("figure2", "C_primary_labels.tsv"),
    output:
        heatmap=result("figure2", "D_effect_heatmap.tsv"),
        concordance=result("figure2", "E_concordance.tsv"),
        ranking=result("figure2", "intermediate/gene_ranking.tsv"),
    params:
        script=script("figure2", "gene_panel_effects.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" analysis="{input.analysis}" '
        'secondary="{input.secondary}" run_metadata="{input.run_metadata}" '
        'nested_predictions="{input.nested_predictions}" sv_predictions="{input.sv_predictions}" '
        'sv_labels="{input.sv_labels}" sv_panels="{input.sv_panels}" labels="{input.labels}" '
        'heatmap="{output.heatmap}" concordance="{output.concordance}" ranking="{output.ranking}"'


rule figure2_primary_label_classifier:
    input:
        ok=FIG2_OK,
        analysis=data("shared/gut_reference_analysis.rds"),
        secondary=data("figure2/classifier_inputs/secondary_validation.rds"),
        core=data(FIG2_RUN + "candidate_panel_core_result.rds"),
        ranking=result("figure2", "intermediate/gene_ranking.tsv"),
        nested_panels=data(FIG2_RUN + "nested_panels.tsv.gz"),
        nested_predictions=data(FIG2_RUN + "nested_predictions.tsv.gz"),
        sv_predictions=data(FIG2_RUN + "sv_predictions.tsv.gz"),
        sv_labels=data(FIG2_RUN + "sv_labels.tsv.gz"),
    output:
        calls=result("figure2", "primary_only_training_calls.tsv"),
        scores=result("figure2", "intermediate/primary_only_cell_scores.tsv"),
    params:
        script=script("figure2", "primary_label_classifier.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" analysis="{input.analysis}" '
        'secondary="{input.secondary}" core="{input.core}" ranking="{input.ranking}" '
        'nested_panels="{input.nested_panels}" nested_predictions="{input.nested_predictions}" '
        'sv_predictions="{input.sv_predictions}" sv_labels="{input.sv_labels}" calls="{output.calls}" '
        'scores="{output.scores}"'


rule figure2_hio_grafts:
    input:
        ok=FIG2_OK,
        counts=data("shared/gse214852_hio/counts.mtx"),
        genes=data("shared/gse214852_hio/genes.txt"),
        cells=data("shared/gse214852_hio/cells.tsv"),
        hio=result("figure2", "intermediate/hio_four_marker_calls.tsv"),
        analysis=data("shared/gut_reference_analysis.rds"),
        ranking=result("figure2", "intermediate/gene_ranking.tsv"),
    output:
        bc=result("figure2", "hio_BC_classifier_calls_by_label.tsv"),
        d=result("figure2", "hio_D_hio_effects.tsv"),
        e=result("figure2", "hio_E_hio_summary.tsv"),
        scores=result("figure2", "intermediate/hio_cell_scores.tsv"),
    params:
        script=script("figure2", "hio_grafts.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" counts="{input.counts}" genes="{input.genes}" '
        'cells="{input.cells}" hio="{input.hio}" analysis="{input.analysis}" ranking="{input.ranking}" '
        'bc="{output.bc}" d="{output.d}" e="{output.e}" scores="{output.scores}"'


rule figure2_c_marker_detection:
    input:
        ok=FIG2_OK,
        cells=result("figure2", "intermediate/heldout_cells.tsv"),
        analysis=data("shared/gut_reference_analysis.rds"),
        secondary=data("figure2/classifier_inputs/secondary_validation.rds"),
        counts=data("shared/gse214852_hio/counts.mtx"),
        genes=data("shared/gse214852_hio/genes.txt"),
        hio_cells=data("shared/gse214852_hio/cells.tsv"),
        hio=result("figure2", "intermediate/hio_four_marker_calls.tsv"),
        labels=result("figure2", "C_primary_labels.tsv"),
        bc=result("figure2", "hio_BC_classifier_calls_by_label.tsv"),
    output:
        result("figure2", "C_marker_detection.tsv"),
    params:
        script=script("figure2", "c_marker_detection.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" analysis="{input.analysis}" '
        'secondary="{input.secondary}" counts="{input.counts}" genes="{input.genes}" '
        'hio_cells="{input.hio_cells}" hio="{input.hio}" labels="{input.labels}" bc="{input.bc}" out="{output}"'


rule figure2_d_detection:
    input:
        ok=FIG2_OK,
        analysis=data("shared/gut_reference_analysis.rds"),
        secondary=data("figure2/classifier_inputs/secondary_validation.rds"),
        sv_labels=data(FIG2_RUN + "sv_labels.tsv.gz"),
        counts=data("shared/gse214852_hio/counts.mtx"),
        genes=data("shared/gse214852_hio/genes.txt"),
        cells=data("shared/gse214852_hio/cells.tsv"),
        hio=result("figure2", "intermediate/hio_four_marker_calls.tsv"),
        heatmap=result("figure2", "D_effect_heatmap.tsv"),
        hio_d=result("figure2", "hio_D_hio_effects.tsv"),
    output:
        result("figure2", "D_detection.tsv"),
    params:
        script=script("figure2", "d_detection.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" analysis="{input.analysis}" secondary="{input.secondary}" '
        'sv_labels="{input.sv_labels}" counts="{input.counts}" genes="{input.genes}" cells="{input.cells}" '
        'hio="{input.hio}" heatmap="{input.heatmap}" hio_d="{input.hio_d}" out="{output}"'


rule figure2_e_cells:
    input:
        ok=FIG2_OK,
        cells=result("figure2", "intermediate/heldout_cells.tsv"),
        scores=result("figure2", "intermediate/primary_only_cell_scores.tsv"),
        hio=result("figure2", "intermediate/hio_four_marker_calls.tsv"),
        hio_scores=result("figure2", "intermediate/hio_cell_scores.tsv"),
        hio_cells=data("shared/gse214852_hio/cells.tsv"),
    output:
        result("figure2", "E_cells.tsv"),
    params:
        script=script("figure2", "e_cells.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" scores="{input.scores}" hio="{input.hio}" '
        'hio_scores="{input.hio_scores}" hio_cells="{input.hio_cells}" out="{output}"'
