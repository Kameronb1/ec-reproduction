# Figure 3 supplement (figure3_supp): rules that recreate the tables behind each panel.
# Rules marked resources heavy_r=1 load the 2.2 GB SC-islet count matrix or the intestinal cell matrix; with
# --resources heavy_r=1 they run one at a time, so only one of them is in memory at a time.
# Script and folder prefixes s1_ to s9_ only group files; they are not figure numbers.


S1_COUNTS = dict(
    meta=data("shared/sc_endocrine_counts/metadata.tsv"),
    genes=data("shared/sc_endocrine_counts/genes.txt"),
    indptr=data("shared/sc_endocrine_counts/counts_indptr.i32"),
    indices=data("shared/sc_endocrine_counts/counts_indices.i32"),
    counts=data("shared/sc_endocrine_counts/counts_data.i32"),
)
S1_ARGS = 'lib="{input.lib}" meta="{input.meta}" genes="{input.genes}" indptr="{input.indptr}" indices="{input.indices}" data="{input.counts}"'


rule f3s_s1_marker_limma:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s1_marker_limma.R"),
        lib=script("figure3_supp", "s1_lib.R"),
        **S1_COUNTS,
    output:
        top20=result("figure3_supp", "s1_marker_limma/top20.tsv"),
        effects=result("figure3_supp", "s1_marker_limma/heatmap_effects.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{input.script}" ' + S1_ARGS + ' top20="{output.top20}" effects="{output.effects}"'


rule f3s_s1_marker_stability:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s1_marker_stability.R"),
        lib=script("figure3_supp", "s1_lib.R"),
        **S1_COUNTS,
    output:
        a=result("figure3_supp", "s1_marker_stability/a_fold_ranks.tsv.gz"),
        c=result("figure3_supp", "s1_marker_stability/c_matched_random.tsv.gz"),
        swaps_input=result("figure3_supp", "s1_marker_stability/swaps_input.rds"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{input.script}" ' + S1_ARGS + ' a="{output.a}" c="{output.c}" swaps_input="{output.swaps_input}"'


# About 1,000 small limma fits on the sample-by-class means; light on memory, slow on one core (workers = R processes).
rule f3s_s1_label_swaps:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s1_label_swaps.R"),
        lib=script("figure3_supp", "s1_lib.R"),
        swaps_input=result("figure3_supp", "s1_marker_stability/swaps_input.rds"),
    output:
        b=result("figure3_supp", "s1_marker_stability/b_label_swaps.tsv.gz"),
    threads: 2
    shell:
        '"{RSCRIPT}" "{input.script}" lib="{input.lib}" swaps_input="{input.swaps_input}" workers={threads} b="{output.b}"'


# Figure S6: classifier gene list x training set swap grid (cells classified as EC per SC-islet test study, A)
# and every SC-islet endocrine cell's four-marker and 20-gene scores and classifications (B).
# The TPH1-only classifier's fit comes from the Figure 2 supplement's rules (loaded with this figure's, see Snakefile).
rule f3s_s2_classifier_swap:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s2_classifier_swap.R"),
        lib=script("figure3_supp", "s2_lib.R"),
        gut_cells=data("figure3_supp/s2/gut_training_cells.tsv.gz"),
        gut_expr=data("figure3_supp/s2/gut_training_panel_gene_expression.tsv.gz"),
        sc_meta=data("shared/sc_endocrine_counts/metadata.tsv"),
        sc_counts=data("figure3_supp/s2/sc_panel_gene_counts.tsv.gz"),
        sc_graphs=data("figure3_supp/s2/sc_graphs.rds"),
        markers=data("shared/four_marker_classifier/fitted_constants.json"),
        tph1_fit=result("figure2_supp", "tph1_only/four_marker_fit.json"),
        gut_genes=data("shared/figure2_panel_d/D_effect_heatmap.tsv"),
        sc_genes=result("figure3_supp", "s1_marker_limma/top20.tsv"),
    output:
        calls=result("figure3_supp", "s2_classifier_swap/calls_by_study.tsv"),
        cells=result("figure3_supp", "s2_classifier_swap/sc_cells.tsv.gz"),
    shell:
        '"{RSCRIPT}" "{input.script}" lib="{input.lib}" gut_cells="{input.gut_cells}" gut_expr="{input.gut_expr}" '
        'sc_meta="{input.sc_meta}" sc_counts="{input.sc_counts}" sc_graphs="{input.sc_graphs}" markers="{input.markers}" '
        'tph1_fit="{input.tph1_fit}" gut_genes="{input.gut_genes}" sc_genes="{input.sc_genes}" out="{output.calls}" cells="{output.cells}"'


# ---- s3_ scripts: Figure S7A (per-donor identity reference), S7B (signature concordance), S9 (artifact sensitivity) ----

S3_PB = "shared/ec_identity_pseudobulk/"
S3_PB_FILES = [S3_PB + t + s for t in ("gut_primary", "sc", "induced", "hio_mccauley_endocrine1")
               for s in ("_sums.tsv.gz", "_detected.tsv.gz", "_groups.tsv")]
S3_SC_COUNTS = ["shared/sc_endocrine_counts/" + f for f in
                ("metadata.tsv", "genes.txt", "counts_indptr.i32", "counts_indices.i32", "counts_data.i32")]
_S3_INPUTS = [r["file"] for r in read_tsv(REPO / "inputs" / "figure3_supp.tsv")]
S3_SC_SAMPLES = [f for f in _S3_INPUTS if f.startswith("shared/sc_sic_sample_pseudobulk/")]
S3_GSE199636 = [f for f in _S3_INPUTS if f.startswith("shared/gse199636_sic_pseudobulk/")]
S3_SOUPX_CAPTURES = [f for f in _S3_INPUTS if f.startswith("figure3_supp/s3/gut_soupx_captures/")]


rule f3s_s3_panel_c_per_donor:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s3_03_panel_c_per_donor.py"),
        lib=script("figure3_supp", "s3_lib.py"),
        pb=[data(f) for f in S3_PB_FILES],
        sc=[data(f) for f in S3_SC_COUNTS],
    output:
        dots=result("figure3_supp", "s3_per_donor/dots_pooled_vs_balanced.tsv"),
    params:
        pb=data("shared/ec_identity_pseudobulk"),
        sc=data("shared/sc_endocrine_counts"),
    resources:
        heavy_r=1,
    shell:
        '"{PYTHON}" "{input.script}" lib="{input.lib}" pb="{params.pb}" sc_counts="{params.sc}" out="{output.dots}"'


rule f3s_s3_signature_scatter:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s3_06_signature_scatter_five.py"),
        lib=script("figure3_supp", "s3_lib.py"),
        pb=[data(S3_PB + "hio_mccauley_endocrine1" + s) for s in ("_sums.tsv.gz", "_detected.tsv.gz", "_groups.tsv")],
        gut=data("shared/ec_identity/gut_source_fold_changes.tsv.gz"),
        signatures=data("shared/ec_identity/q1_gene_signatures.tsv.gz"),
    output:
        genes=result("figure3_supp", "s3_signature/gene_signatures.tsv.gz"),
        concordance=result("figure3_supp", "s3_signature/concordance.tsv"),
    params:
        pb=data("shared/ec_identity_pseudobulk"),
    shell:
        '"{PYTHON}" "{input.script}" lib="{input.lib}" pb="{params.pb}" gut="{input.gut}" signatures="{input.signatures}" '
        'out_genes="{output.genes}" out_concordance="{output.concordance}"'


rule f3s_s3_share_two_fold:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s3_09_share_two_fold.py"),
        genes=result("figure3_supp", "s3_signature/gene_signatures.tsv.gz"),
    output:
        share=result("figure3_supp", "s3_signature/share_two_fold.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" genes="{input.genes}" out="{output.share}"'


rule f3s_s3_fig3d_artifact_checks:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s3_16_fig3d_artifact_checks.R"),
        analysis=data("shared/gut_reference_analysis.rds"),
        sc_files=[data(f) for f in S3_SC_SAMPLES],
        g6_files=[data(f) for f in S3_GSE199636],
        fig3d=data("shared/sc_sic_vs_gut_donor_block/contrast_donor_block_levels.tsv.gz"),
        captures=[data(f) for f in S3_SOUPX_CAPTURES],
        sc_soupx=data("figure3_supp/s3/sic_counts_raw_and_soupx.rds"),
        gut_doublets=data("figure3_supp/s3/cells.tsv.gz"),
        sc_doublets=data("figure3_supp/s3/sc_sic_cells.tsv.gz"),
    output:
        genes=result("figure3_supp", "s3_artifact_checks/genes_before_after.tsv.gz"),
    params:
        sc_groups=data("shared/sc_sic_sample_pseudobulk/groups.tsv.gz"),
        sc_pb=data("shared/sc_sic_sample_pseudobulk"),
        g6_groups=data("shared/gse199636_sic_pseudobulk/groups.tsv"),
        g6_pb=data("shared/gse199636_sic_pseudobulk"),
        captures=data("figure3_supp/s3/gut_soupx_captures"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{input.script}" analysis="{input.analysis}" sc_groups="{params.sc_groups}" sc_pb="{params.sc_pb}" '
        'gse199636_groups="{params.g6_groups}" gse199636_pb="{params.g6_pb}" fig3d="{input.fig3d}" '
        'gut_soupx="{params.captures}" sc_soupx="{input.sc_soupx}" gut_doublets="{input.gut_doublets}" '
        'sc_doublets="{input.sc_doublets}" out="{output.genes}"'


# Figure S8: HIO graft EC versus SC-SIC versus primary intestinal EC (A) and the enrichment of its gene patterns (B).

rule f3s_s4_three_way_patterns:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s4_hio_three_way_patterns.R"),
        gut_sums=data("shared/ec_identity_pseudobulk/gut_primary_sums.tsv.gz"),
        gut_groups=data("shared/ec_identity_pseudobulk/gut_primary_groups.tsv"),
        sc_sums=data("shared/ec_identity_pseudobulk/sc_sums.tsv.gz"),
        sc_groups=data("shared/ec_identity_pseudobulk/sc_groups.tsv"),
        hio_sums=data("shared/ec_identity_pseudobulk/hio_mccauley_endocrine1_sums.tsv.gz"),
        hio_groups=data("shared/ec_identity_pseudobulk/hio_mccauley_endocrine1_groups.tsv"),
        roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
    output:
        out=result("figure3_supp", "s4_hio/genes_all.tsv.gz"),
    shell:
        '"{RSCRIPT}" "{input.script}" gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" sc_sums="{input.sc_sums}" '
        'sc_groups="{input.sc_groups}" hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" roster="{input.roster}" '
        'out="{output.out}"'


rule f3s_s4_enrichment:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s4_hio_enrichment.R"),
        genes_all=result("figure3_supp", "s4_hio/genes_all.tsv.gz"),
        msigdb_zip=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all=result("figure3_supp", "s4_hio/enrichment/enrichment_all_terms.tsv"),
        called=result("figure3_supp", "s4_hio/enrichment/enrichment_called.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" genes_all="{input.genes_all}" msigdb_zip="{input.msigdb_zip}" '
        'all="{output.all}" called="{output.called}"'


# MSigDB C8 cell-type enrichment of the six HIO pattern gene lists (Figure S8B).
rule f3s_s4_celltype_signatures:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s4_hio_celltype_signatures.R"),
        genes_all=result("figure3_supp", "s4_hio/genes_all.tsv.gz"),
        msigdb_zip=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all=result("figure3_supp", "s4_hio/enrichment_celltype/celltype_all_terms.tsv"),
        called=result("figure3_supp", "s4_hio/enrichment_celltype/celltype_called.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" genes_all="{input.genes_all}" msigdb_zip="{input.msigdb_zip}" '
        'all="{output.all}" called="{output.called}"'


rule f3s_s6_gene_group_tables:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s6_gene_group_tables.R"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        gut_sums=data("shared/ec_identity_pseudobulk/gut_primary_sums.tsv.gz"),
        gut_groups=data("shared/ec_identity_pseudobulk/gut_primary_groups.tsv"),
        gut_roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
        sc_sums=data("shared/ec_identity_pseudobulk/sc_sums.tsv.gz"),
        sc_groups=data("shared/ec_identity_pseudobulk/sc_groups.tsv"),
    output:
        s6=result("figure3_supp", "s6_gene_groups/g1_g2_heatmap.tsv"),
        s7=result("figure3_supp", "s6_gene_groups/g3_heatmap.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" groups="{input.groups}" gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" gut_roster="{input.gut_roster}" '
        'sc_sums="{input.sc_sums}" sc_groups="{input.sc_groups}" s6="{output.s6}" s7="{output.s7}"'


# Figure S13: Figure 3 gene groups with and without their mixed genes; pathways (A) and C8 cell types (B).

rule f3s_s8_mixed_gene_enrichment:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s8_mixed_gene_enrichment.R"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        mixed=data("shared/gene_groups/mixed_genes.tsv"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        msigdb_zip=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all=result("figure3_supp", "s8_mixed/pathways/enrichment_all_terms.tsv.gz"),
        called=result("figure3_supp", "s8_mixed/pathways/enrichment_called.tsv"),
        compare=result("figure3_supp", "s8_mixed/pathways/group_vs_group_with_mixed.tsv"),
        lists=result("figure3_supp", "s8_mixed/pathways/lists.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" groups="{input.groups}" mixed="{input.mixed}" universe="{input.universe}" '
        'msigdb_zip="{input.msigdb_zip}" all="{output.all}" called="{output.called}" compare="{output.compare}" lists="{output.lists}"'


# MSigDB C8 cell-type enrichment of each gene group alone and with its mixed genes (Figure S13B).
rule f3s_s8_mixed_gene_celltypes:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s8_mixed_gene_celltypes.R"),
        lists=result("figure3_supp", "s8_mixed/pathways/lists.tsv"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        msigdb_zip=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all=result("figure3_supp", "s8_mixed/celltypes/celltype_all_terms.tsv"),
        called=result("figure3_supp", "s8_mixed/celltypes/celltype_called.tsv"),
        compare=result("figure3_supp", "s8_mixed/celltypes/group_vs_group_with_mixed.tsv"),
        by_class=result("figure3_supp", "s8_mixed/celltypes/called_terms_by_class.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" lists="{input.lists}" universe="{input.universe}" msigdb_zip="{input.msigdb_zip}" '
        'all="{output.all}" called="{output.called}" compare="{output.compare}" by_class="{output.by_class}"'


PB_S9 = "shared/ec_identity_pseudobulk/"


rule f3s_s9_brain_similarity:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s9_brain_similarity.R"),
        gut_sums=data(PB_S9 + "gut_primary_sums.tsv.gz"),
        gut_groups=data(PB_S9 + "gut_primary_groups.tsv"),
        sc_sums=data(PB_S9 + "sc_sums.tsv.gz"),
        sc_groups=data(PB_S9 + "sc_groups.tsv"),
        hio_sums=data(PB_S9 + "hio_mccauley_endocrine1_sums.tsv.gz"),
        hio_groups=data(PB_S9 + "hio_mccauley_endocrine1_groups.tsv"),
        ind_sums=data(PB_S9 + "induced_sums.tsv.gz"),
        ind_groups=data(PB_S9 + "induced_groups.tsv"),
        islet_sums=data(PB_S9 + "islet_sums.tsv.gz"),
        islet_groups=data(PB_S9 + "islet_groups.tsv"),
        gut_roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
        brain_sums=data("figure3_supp/s9/siletti2023_neurons_sums.tsv.gz"),
        brain_groups=data("figure3_supp/s9/siletti2023_neurons_groups.tsv"),
        islet_roster=data("shared/adult_islet_roster/roster.json"),
    output:
        similarity=result("figure3_supp", "s9_brain/similarity_units.tsv"),
        units=result("figure3_supp", "s9_brain/raphe_units.tsv"),
        summary=result("figure3_supp", "s9_brain/summary_by_group.tsv"),
        ser_da_genes=result("figure3_supp", "s9_brain/ser_da_genes.tsv"),
        ser_genes=result("figure3_supp", "s9_brain/ser_genes_in_sc_sic.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" sc_sums="{input.sc_sums}" '
        'sc_groups="{input.sc_groups}" hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" ind_sums="{input.ind_sums}" '
        'ind_groups="{input.ind_groups}" islet_sums="{input.islet_sums}" islet_groups="{input.islet_groups}" '
        'gut_roster="{input.gut_roster}" brain_sums="{input.brain_sums}" brain_groups="{input.brain_groups}" islet_roster="{input.islet_roster}" '
        'similarity="{output.similarity}" units="{output.units}" summary="{output.summary}" '
        'ser_da_genes="{output.ser_da_genes}" ser_genes="{output.ser_genes}"'


rule f3s_s9_stac_by_sample:
    input:
        ok=str(RESULTS / "figure3_supp_inputs_verified.txt"),
        script=script("figure3_supp", "s9_stac_by_sample.R"),
        gut_sums=data(PB_S9 + "gut_primary_sums.tsv.gz"),
        gut_detected=data(PB_S9 + "gut_primary_detected.tsv.gz"),
        gut_groups=data(PB_S9 + "gut_primary_groups.tsv"),
        hio_sums=data(PB_S9 + "hio_mccauley_endocrine1_sums.tsv.gz"),
        hio_detected=data(PB_S9 + "hio_mccauley_endocrine1_detected.tsv.gz"),
        hio_groups=data(PB_S9 + "hio_mccauley_endocrine1_groups.tsv"),
        ind_sums=data(PB_S9 + "induced_sums.tsv.gz"),
        ind_detected=data(PB_S9 + "induced_detected.tsv.gz"),
        ind_groups=data(PB_S9 + "induced_groups.tsv"),
        sc_sums=data(PB_S9 + "sc_sums.tsv.gz"),
        sc_detected=data(PB_S9 + "sc_detected.tsv.gz"),
        sc_groups=data(PB_S9 + "sc_groups.tsv"),
        islet_sums=data(PB_S9 + "islet_sums.tsv.gz"),
        islet_detected=data(PB_S9 + "islet_detected.tsv.gz"),
        islet_groups=data(PB_S9 + "islet_groups.tsv"),
        gut_roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
        siletti_sums=data("figure3_supp/s9/siletti2023_neurons_sums.tsv.gz"),
        siletti_detected=data("figure3_supp/s9/siletti2023_neurons_detected.tsv.gz"),
        siletti_groups=data("figure3_supp/s9/siletti2023_neurons_groups.tsv"),
        braun_sums=data("figure3_supp/s9/braun2023_neurons_sums.tsv.gz"),
        braun_detected=data("figure3_supp/s9/braun2023_neurons_detected.tsv.gz"),
        braun_groups=data("figure3_supp/s9/braun2023_neurons_groups.tsv"),
        islet_roster=data("shared/adult_islet_roster/roster.json"),
    output:
        by_sample=result("figure3_supp", "s9_stac/stac_by_sample.tsv"),
        drawn=result("figure3_supp", "s9_stac/stac_by_sample_drawn.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" gut_sums="{input.gut_sums}" gut_detected="{input.gut_detected}" gut_groups="{input.gut_groups}" '
        'hio_sums="{input.hio_sums}" hio_detected="{input.hio_detected}" hio_groups="{input.hio_groups}" '
        'ind_sums="{input.ind_sums}" ind_detected="{input.ind_detected}" ind_groups="{input.ind_groups}" '
        'sc_sums="{input.sc_sums}" sc_detected="{input.sc_detected}" sc_groups="{input.sc_groups}" '
        'islet_sums="{input.islet_sums}" islet_detected="{input.islet_detected}" islet_groups="{input.islet_groups}" '
        'gut_roster="{input.gut_roster}" siletti_sums="{input.siletti_sums}" siletti_detected="{input.siletti_detected}" '
        'siletti_groups="{input.siletti_groups}" braun_sums="{input.braun_sums}" braun_detected="{input.braun_detected}" '
        'braun_groups="{input.braun_groups}" islet_roster="{input.islet_roster}" '
        'by_sample="{output.by_sample}" drawn="{output.drawn}"'
