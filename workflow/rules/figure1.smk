# Figure 1: study sources and ages (A), study roles (B), per-study display embeddings (C),
# serotonin transcript co-detection and the number of the three genes per cell (D), islet hormone
# transcripts in the cells with TPH1 detected (E).

FIG1_OK = str(RESULTS / "figure1_inputs_verified.txt")


rule figure1_native_age_chips:
    input:
        ok=FIG1_OK,
        ages=data("figure1/native_age_records.tsv"),
    output:
        result("figure1", "native_age_chips.tsv"),
    params:
        script=script("figure1", "native_age_chips.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" ages="{input.ages}" out="{output}"'


rule figure1_study_roles:
    input:
        ok=FIG1_OK,
        roles=data("figure1/study_roles.tsv"),
    output:
        result("figure1", "study_roles_drawn.tsv"),
    params:
        script=script("figure1", "study_roles.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" roles="{input.roles}" out="{output}"'


rule figure1_display_embedding:
    input:
        ok=FIG1_OK,
        cells=data("figure1/endocrine_display_cells.tsv.gz"),
        calls=data("figure1/fig2_four_marker_calls.tsv.gz"),
        roles=result("figure1", "study_roles_drawn.tsv"),
    output:
        result("figure1", "display_embedding_cells.tsv.gz"),
    params:
        script=script("figure1", "display_embedding.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" calls="{input.calls}" '
        'roles="{input.roles}" out="{output}"'


rule figure1_serotonin_codetection:
    input:
        ok=FIG1_OK,
        cells=data("figure1/endocrine_display_cells.tsv.gz"),
        sc=data("figure1/sc_islet_soupx_serotonin_genes.tsv.gz"),
        gut_ec=data("figure1/gut_ec_soupx_serotonin_genes.tsv.gz"),
        gut_other=data("figure1/gut_other_eec_soupx_serotonin_genes.tsv.gz"),
        adult=data("figure1/adult_islet_soupx_serotonin_genes.tsv.gz"),
    output:
        result("figure1", "serotonin_codetection.tsv"),
    params:
        script=script("figure1", "serotonin_codetection.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" sc="{input.sc}" '
        'gut_ec="{input.gut_ec}" gut_other="{input.gut_other}" adult="{input.adult}" out="{output}"'


rule figure1_serotonin_gene_count:
    input:
        ok=FIG1_OK,
        cells=data("figure1/endocrine_display_cells.tsv.gz"),
        sc=data("figure1/sc_islet_soupx_serotonin_genes.tsv.gz"),
        gut_ec=data("figure1/gut_ec_soupx_serotonin_genes.tsv.gz"),
        gut_other=data("figure1/gut_other_eec_soupx_serotonin_genes.tsv.gz"),
        adult=data("figure1/adult_islet_soupx_serotonin_genes.tsv.gz"),
        codetection=result("figure1", "serotonin_codetection.tsv"),
    output:
        gene_count=result("figure1", "serotonin_gene_count.tsv"),
        hormones=result("figure1", "hormones_in_tph1_cells.tsv"),
    params:
        script=script("figure1", "serotonin_gene_count.R"),
    shell:
        '"{RSCRIPT}" "{params.script}" cells="{input.cells}" sc="{input.sc}" '
        'gut_ec="{input.gut_ec}" gut_other="{input.gut_other}" adult="{input.adult}" codetection="{input.codetection}" '
        'out_count="{output.gene_count}" out_hormones="{output.hormones}"'
