# Figure 1: the six SC-islet studies on their display embeddings and the SC-SIC share of each (C), serotonin transcript
# co-detection and the number of the three genes per cell (D), and the four-marker score of SC-islet cells against the
# classifier's cutoff (E). The islet hormone table drawn in Supplementary Figure S1A is written here too, by the script
# that writes D's gene counts. E reads the classifier results of the Supplementary Figure S2-S3 unit.

FIG1_OK = str(RESULTS / "figure1_inputs_verified.txt")
FIG1_CELLS = data("figure1/endocrine_display_cells.tsv.gz")
FIG1_SOUPX = dict(
    sc=data("figure1/sc_islet_soupx_serotonin_genes.tsv.gz"),
    gut_ec=data("figure1/gut_ec_soupx_serotonin_genes.tsv.gz"),
    gut_other=data("figure1/gut_other_eec_soupx_serotonin_genes.tsv.gz"),
    adult=data("figure1/adult_islet_soupx_serotonin_genes.tsv.gz"),
)
FIG1_SOUPX_ARGS = 'sc="{input.sc}" gut_ec="{input.gut_ec}" gut_other="{input.gut_other}" adult="{input.adult}"'


rule figure1_display_embedding:
    input:
        ok=FIG1_OK,
        cells=FIG1_CELLS,
    output:
        cells=result("figure1", "display_embedding_cells.tsv.gz"),
        share=result("figure1", "sc_sic_share.tsv"),
    params:
        script=script("figure1", "display_embedding.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" out="{output.cells}" share="{output.share}"'


rule figure1_serotonin_codetection:
    input:
        ok=FIG1_OK,
        cells=FIG1_CELLS,
        **FIG1_SOUPX,
    output:
        result("figure1", "serotonin_codetection.tsv"),
    params:
        script=script("figure1", "serotonin_codetection.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" ' + FIG1_SOUPX_ARGS + ' out="{output}"'


rule figure1_serotonin_gene_count:
    input:
        ok=FIG1_OK,
        cells=FIG1_CELLS,
        **FIG1_SOUPX,
        codetection=result("figure1", "serotonin_codetection.tsv"),
    output:
        gene_count=result("figure1", "serotonin_gene_count.tsv"),
        groups=result("figure1", "hormones_by_group.tsv"),
    params:
        script=script("figure1", "serotonin_gene_count.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" ' + FIG1_SOUPX_ARGS + ' codetection="{input.codetection}" '
        'out_count="{output.gene_count}" out_groups="{output.groups}"'


rule figure1_sc_islet_above_cutoff:
    input:
        ok=FIG1_OK,
        scores=result("figure_s2_s3", "sc_islet_scores_by_cell.tsv.gz"),
        cutoff=result("figure_s2_s3", "sc_islet_cutoff.tsv"),
        constants=result("figure_s2_s3", "fitted_constants.json"),
    output:
        result("figure1", "sc_islet_above_cutoff.tsv"),
    params:
        script=script("figure1", "sc_islet_above_cutoff.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" scores="{input.scores}" cutoff="{input.cutoff}" constants="{input.constants}" '
        'out="{output}"'
