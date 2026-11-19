# Figures S2 and S3: the label-free four-marker EC classifier. Its fit on the 1,164 fitting cells (the seven pediatric
# Crohn's disease donors of Elmentaite 2020 left out; S2A), the cells it classifies as EC within the published labels
# of the four test studies (S2B), the marker-free map and clusters of the fitting cells (S2C-D), classification across
# cutoffs with the leave-one-fitting-group-out refits (S3A), read thinning (S3B), and SC-SICs and other SC endocrine
# cells classified as EC (S3C).
# Three results are read by other figures: fitted_constants.json (the classifier's constants),
# sc_islet_scores_by_cell.tsv.gz and sc_islet_cutoff.tsv (the four-marker score of every SC-islet cell, and the cutoff).
# Rules marked resources heavy_r=1 hold the fitting-cell matrix or a test study's full count matrix in memory; with
# --resources heavy_r=1 they run one at a time.

S2S3_OK = str(RESULTS / "figure_s2_s3_inputs_verified.txt")
S2S3_MATRIX = dict(
    matrix=data("shared/fitting_matrix/fitting_matrix.f64le"),
    features=data("shared/fitting_matrix/features.tsv"),
    cells=data("shared/fitting_matrix/cells.tsv"),
    meta=data("shared/fitting_matrix/cell_metadata.tsv.gz"),
    roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
)
S2S3_MATRIX_ARGS = ('matrix="{input.matrix}" features="{input.features}" cells="{input.cells}" meta="{input.meta}" '
                    'roster="{input.roster}" ')
S2S3_H5AD = [data(f"figure2/heldout_studies/{s}.h5ad") for s in ("gut_cell_atlas", "hickey", "egozi")]
S2S3_H5AD_DIR = str(Path(S2S3_H5AD[0]).parent)
S2S3_HIO_CELLS = data("shared/gse214852_hio/cells.tsv")


def s2s3(name):
    return result("figure_s2_s3", name)


rule figure_s2_s3_marker_fit:
    input:
        ok=S2S3_OK,
        **S2S3_MATRIX,
    output:
        fit=s2s3("fitted_constants.json"),
        cells_out=s2s3("fitting_cell_classes.tsv.gz"),
        counts=s2s3("fitting_class_counts.tsv"),
        hio=s2s3("intermediate/hpsc_cell_classes.tsv.gz"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "marker_fit.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" ' + S2S3_MATRIX_ARGS +
        'fit_out="{output.fit}" cells_out="{output.cells_out}" counts_out="{output.counts}" hio_cells_out="{output.hio}"'


rule figure_s2_s3_four_marker_fit:
    input:
        ok=S2S3_OK,
        constants=s2s3("fitted_constants.json"),
        **S2S3_MATRIX,
    output:
        hist=s2s3("A_histogram.tsv"),
        curves=s2s3("A_mixture_curves.tsv"),
        hio=s2s3("intermediate/hio_four_marker_calls.tsv"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "four_marker_fit.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" ' + S2S3_MATRIX_ARGS +
        'constants="{input.constants}" hist="{output.hist}" curves="{output.curves}" hio="{output.hio}"'


rule figure_s2_s3_cutoff_stability:
    input:
        ok=S2S3_OK,
        **S2S3_MATRIX,
    output:
        s2s3("cutoff_stability_summary.tsv"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "cutoff_stability.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" ' + S2S3_MATRIX_ARGS + 'out="{output}"'


rule figure_s2_s3_umap:
    input:
        ok=S2S3_OK,
        classes=s2s3("fitting_cell_classes.tsv.gz"),
        **S2S3_MATRIX,
    output:
        coords=s2s3("umap_coordinates.tsv.gz"),
        embeddings=s2s3("harmony_embeddings.tsv.gz"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "umap.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" ' + S2S3_MATRIX_ARGS +
        'classes="{input.classes}" coords="{output.coords}" embeddings="{output.embeddings}"'


rule figure_s2_s3_clusters:
    input:
        ok=S2S3_OK,
        classes=s2s3("fitting_cell_classes.tsv.gz"),
        embeddings=s2s3("harmony_embeddings.tsv.gz"),
        **S2S3_MATRIX,
    output:
        summary=s2s3("cluster_summary.tsv"),
        assign=s2s3("cell_clusters.tsv.gz"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "clusters.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" ' + S2S3_MATRIX_ARGS +
        'classes="{input.classes}" embeddings="{input.embeddings}" summary="{output.summary}" assign="{output.assign}"'


rule figure_s2_s3_heldout_marker_calls:
    input:
        ok=S2S3_OK,
        fit=s2s3("fitted_constants.json"),
        h5ad=S2S3_H5AD,
    output:
        labels=s2s3("C_primary_labels.tsv"),
        cells=s2s3("intermediate/heldout_cells.tsv"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "heldout_marker_calls.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" fit="{input.fit}" h5ad_gut_cell_atlas="{input.h5ad[0]}" '
        'h5ad_hickey="{input.h5ad[1]}" h5ad_egozi="{input.h5ad[2]}" labels="{output.labels}" cells="{output.cells}"'


rule figure_s2_s3_hio_calls_by_label:
    input:
        ok=S2S3_OK,
        cells=S2S3_HIO_CELLS,
        hio=s2s3("intermediate/hio_four_marker_calls.tsv"),
    output:
        s2s3("hio_BC_classifier_calls_by_label.tsv"),
    params:
        script=script("figure_s2_s3", "hio_calls_by_label.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" hio="{input.hio}" bc="{output}"'


rule figure_s2_s3_marker_detection:
    input:
        ok=S2S3_OK,
        cells=s2s3("intermediate/heldout_cells.tsv"),
        h5ad=S2S3_H5AD,
        counts=data("shared/gse214852_hio/counts.mtx"),
        genes=data("shared/gse214852_hio/genes.txt"),
        hio_cells=S2S3_HIO_CELLS,
        hio=s2s3("intermediate/hio_four_marker_calls.tsv"),
        labels=s2s3("C_primary_labels.tsv"),
        bc=s2s3("hio_BC_classifier_calls_by_label.tsv"),
    output:
        s2s3("C_marker_detection.tsv"),
    params:
        script=script("figure_s2_s3", "c_marker_detection.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" h5ad_gut_cell_atlas="{input.h5ad[0]}" '
        'h5ad_hickey="{input.h5ad[1]}" h5ad_egozi="{input.h5ad[2]}" counts="{input.counts}" genes="{input.genes}" '
        'hio_cells="{input.hio_cells}" hio="{input.hio}" labels="{input.labels}" bc="{input.bc}" out="{output}"'


rule figure_s2_s3_heldout:
    input:
        ok=S2S3_OK,
        h5ad=S2S3_H5AD,
        fit=s2s3("fitted_constants.json"),
        stability=s2s3("cutoff_stability_summary.tsv"),
        meta=data("shared/fitting_matrix/cell_metadata.tsv.gz"),
        roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
    output:
        scoring=s2s3("intermediate/heldout_scoring.tsv.gz"),
        sweep=s2s3("intermediate/cutoff_sweep.tsv.gz"),
        refits=s2s3("intermediate/transported_refits.tsv"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "heldout.R"),
        h5ad_dir=S2S3_H5AD_DIR,
    shell:
        '"{RSCRIPT}" "{params.script}" h5ad_dir="{params.h5ad_dir}" fit="{input.fit}" stability="{input.stability}" '
        'meta="{input.meta}" roster="{input.roster}" scoring="{output.scoring}" sweep="{output.sweep}" refits="{output.refits}"'


rule figure_s2_s3_hio_four_marker:
    input:
        ok=S2S3_OK,
        fit=s2s3("fitted_constants.json"),
        stability=s2s3("cutoff_stability_summary.tsv"),
        hio_classes=s2s3("intermediate/hpsc_cell_classes.tsv.gz"),
        hio_cells=S2S3_HIO_CELLS,
        refits=s2s3("intermediate/transported_refits.tsv"),
        **S2S3_MATRIX,
    output:
        sweep=s2s3("intermediate/hio_cutoff_sweep.tsv.gz"),
        refits_out=s2s3("intermediate/hio_transported_refits.tsv"),
        depth=s2s3("hio_four_marker_depth.tsv"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "hio_four_marker.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" ' + S2S3_MATRIX_ARGS +
        'fit="{input.fit}" stability="{input.stability}" hio_classes="{input.hio_classes}" hio_cells="{input.hio_cells}" '
        'refits="{input.refits}" sweep="{output.sweep}" refits_out="{output.refits_out}" depth="{output.depth}"'


rule figure_s2_s3_four_marker_depth:
    input:
        ok=S2S3_OK,
        h5ad=S2S3_H5AD,
        fit=s2s3("fitted_constants.json"),
        scoring=s2s3("intermediate/heldout_scoring.tsv.gz"),
    output:
        s2s3("four_marker_depth.tsv"),
    resources:
        heavy_r=1,
    params:
        script=script("figure_s2_s3", "four_marker_depth.R"),
        h5ad_dir=S2S3_H5AD_DIR,
    shell:
        '"{RSCRIPT}" "{params.script}" h5ad_dir="{params.h5ad_dir}" fit="{input.fit}" scoring="{input.scoring}" out="{output}"'


rule figure_s2_s3_s3a_table:
    input:
        ok=S2S3_OK,
        constants=s2s3("fitted_constants.json"),
        sweep=s2s3("intermediate/cutoff_sweep.tsv.gz"),
        refits=s2s3("intermediate/transported_refits.tsv"),
        hio_sweep=s2s3("intermediate/hio_cutoff_sweep.tsv.gz"),
        hio_refits=s2s3("intermediate/hio_transported_refits.tsv"),
    output:
        s2s3("S3A_classified_across_cutoffs.tsv"),
    params:
        script=script("figure_s2_s3", "s3a_table.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" constants="{input.constants}" sweep="{input.sweep}" refits="{input.refits}" '
        'hio_sweep="{input.hio_sweep}" hio_refits="{input.hio_refits}" out="{output}"'


rule figure_s2_s3_sc_islet_scores:
    input:
        ok=S2S3_OK,
        constants=s2s3("fitted_constants.json"),
        metadata=data("shared/sc_endocrine_counts/metadata.tsv"),
        graphs=data("figure3_supp/s2/sc_graphs.rds"),
        counts=data("figure3_supp/s2/sc_panel_gene_counts.tsv.gz"),
        cells=data("figure3/f3_03_scores_mature_endocrine.tsv.gz"),
    output:
        scores=s2s3("sc_islet_scores_by_cell.tsv.gz"),
        cutoff=s2s3("sc_islet_cutoff.tsv"),
        classified=s2s3("classified_by_study.tsv"),
    params:
        script=script("figure_s2_s3", "score_sc_islet.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" constants="{input.constants}" metadata="{input.metadata}" graphs="{input.graphs}" '
        'counts="{input.counts}" cells="{input.cells}" scores="{output.scores}" cutoff="{output.cutoff}" '
        'classified="{output.classified}"'
