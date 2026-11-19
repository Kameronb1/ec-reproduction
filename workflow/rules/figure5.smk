# Figure 5 (panels A-C): the tables behind each panel, and the regulon and CDX2 knockout tables that Supplementary
# Figure S9 draws. Every rule takes the input check as its first input so nothing runs on unverified data.
# The gene groups G1-G3 are Figure 2's result.

F5_OK = str(RESULTS / "figure5_inputs_verified.txt")
PB = "shared/ec_identity_pseudobulk/"
REG = "figure5/regulons/"
KO = "figure5/cdx2_knockout/"
GENE_GROUPS = result("figure2", "gene_groups.tsv.gz")


rule figure5_regulator_rna:
    """Figure 5A: thinned pseudobulk RNA of 34 candidate regulators per cell type, and SC-SIC minus SC-beta per study."""
    input:
        ok=F5_OK,
        script=script("figure5", "regulator_rna.R"),
        sc_sums=data(PB + "sc_sums.tsv.gz"), sc_groups=data(PB + "sc_groups.tsv"),
        gut_sums=data(PB + "gut_primary_sums.tsv.gz"), gut_groups=data(PB + "gut_primary_groups.tsv"),
        islet_sums=data(PB + "islet_sums.tsv.gz"), islet_groups=data(PB + "islet_groups.tsv"),
        hio_sums=data(PB + "hio_mccauley_endocrine1_sums.tsv.gz"), hio_groups=data(PB + "hio_mccauley_endocrine1_groups.tsv"),
    output:
        by_sample=result("figure5", "regulator_rna/tf_expression_by_sample.tsv"),
        pooled=result("figure5", "regulator_rna/tf_expression_pooled_sparse_gut.tsv"),
        panel=result("figure5", "regulator_rna/f5_5a_regulator_rna.tsv"),
        hio_rows=result("figure5", "regulator_rna/f5_5a_regulator_rna__hio_mccauley_rows.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" sc_sums="{input.sc_sums}" sc_groups="{input.sc_groups}" '
        'gut_primary_sums="{input.gut_sums}" gut_primary_groups="{input.gut_groups}" '
        'islet_sums="{input.islet_sums}" islet_groups="{input.islet_groups}" '
        'hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" '
        'by_sample="{output.by_sample}" pooled="{output.pooled}" panel="{output.panel}" hio_rows="{output.hio_rows}"'


rule figure5_regulon_selection:
    """Figure 5B: candidate regulators whose regulons are reproducible in SC-islets and intestine."""
    input:
        ok=F5_OK,
        script=script("figure5", "regulon_selection.py"),
        regulons=data(REG + "network_regulons.tsv.gz"),
        training=data(REG + "network_training_studies.tsv"),
        effects=data(REG + "regulon_target_effects.tsv.gz"),
        gene_bins=data("figure5/regulator_rna_targets_main/matched_null_gene_bins.tsv"),
        gene_groups=GENE_GROUPS,
        symbol_updates=script("figure5", "target_symbol_updates.tsv"),
    output:
        all_programs=result("figure5", "regulons/all_programs.tsv"),
        factors=result("figure5", "regulons/factors.tsv"),
        program_values=result("figure5", "regulons/program_values.tsv"),
        gut_parts=result("figure5", "regulons/gut_parts.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" regulons="{input.regulons}" training="{input.training}" effects="{input.effects}" '
        'gene_bins="{input.gene_bins}" gene_groups="{input.gene_groups}" symbol_updates="{input.symbol_updates}" '
        'all_programs="{output.all_programs}" factors="{output.factors}" program_values="{output.program_values}" '
        'gut_parts="{output.gut_parts}"'


rule figure5_regulon_target_placement:
    """Figure S9B: where the pooled targets of the selected candidate regulators fall among G2, G3 and beta-side genes."""
    input:
        ok=F5_OK,
        script=script("figure5", "regulon_target_placement.py"),
        factors=result("figure5", "regulons/factors.tsv"),
        all_programs=result("figure5", "regulons/all_programs.tsv"),
        regulons=data(REG + "network_regulons.tsv.gz"),
        gene_bins=data("figure5/regulator_rna_targets_main/matched_null_gene_bins.tsv"),
        gene_groups=GENE_GROUPS,
        sums=data(PB + "sc_sums.tsv.gz"),
        groups=data(PB + "sc_groups.tsv"),
        motif_families=data("figure5/chromatin_motif_screen_g2g3/motif_family_definitions.tsv"),
        motif_hits=data("figure5/chromatin_motif_screen_g2g3/replicated_family_hits.tsv"),
        symbol_updates=script("figure5", "target_symbol_updates.tsv"),
    output:
        tests=result("figure5", "regulons/tests.tsv"),
        factors_display=result("figure5", "regulons/factors_display.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" factors="{input.factors}" all_programs="{input.all_programs}" regulons="{input.regulons}" '
        'gene_bins="{input.gene_bins}" gene_groups="{input.gene_groups}" pseudobulk_sums="{input.sums}" '
        'pseudobulk_groups="{input.groups}" motif_families="{input.motif_families}" motif_hits="{input.motif_hits}" '
        'symbol_updates="{input.symbol_updates}" tests="{output.tests}" factors_display="{output.factors_display}"'


rule figure5_regulon_mccauley_check:
    """Figure S9A, HIO graft EC column: the selected regulons scored in McCauley 2023 HIO graft EC on the published labels."""
    input:
        ok=F5_OK,
        script=script("figure5", "regulon_mccauley_check.py"),
        program_values=result("figure5", "regulons/program_values.tsv"),
        effects=data(REG + "regulon_target_effects.tsv.gz"),
        training=data(REG + "network_training_studies.tsv"),
        cells=data(REG + "mccauley2023_cells.tsv"),
        auc=data(REG + "mccauley2023_regulon_auc.tsv.gz"),
    output:
        out=result("figure5", "regulons/mccauley2023_values.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" program_values="{input.program_values}" effects="{input.effects}" '
        'training="{input.training}" cells="{input.cells}" auc="{input.auc}" out="{output.out}"'


rule figure5_regulon_supplement_values:
    """Figure S9A (intestinal columns) and S9B: intestinal study means per regulon, and G2/G3 among pooled targets."""
    input:
        ok=F5_OK,
        script=script("figure5", "regulon_supplement_values.R"),
        values=result("figure5", "regulons/program_values.tsv"),
        parts=result("figure5", "regulons/gut_parts.tsv"),
        factors=result("figure5", "regulons/factors.tsv"),
        tests=result("figure5", "regulons/tests.tsv"),
        fdisp=result("figure5", "regulons/factors_display.tsv"),
    output:
        out=result("figure5", "regulons/supplement_values.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" values="{input.values}" parts="{input.parts}" factors="{input.factors}" '
        'tests="{input.tests}" fdisp="{input.fdisp}" out="{output.out}"'


rule figure5_cdx2ko_group_test:
    """Figure 5C: knockout minus wild type on the gene groups in SC-SIC, SC-beta and SC-alpha against 500
    expression-matched random gene sets (about 40 s)."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx2ko_group_test.py"),
        counts=data(KO + "recounted_counts.h5ad"),
        labels=data(KO + "transferred_labels.tsv.gz"),
        groups=GENE_GROUPS,
    output:
        effects=result("figure5", "cdx2_knockout/group_test/group_effects.tsv"),
        nuclei=result("figure5", "cdx2_knockout/group_test/nuclei.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" counts="{input.counts}" labels="{input.labels}" groups="{input.groups}" '
        'effects="{output.effects}" nuclei="{output.nuclei}"'


rule figure5_cdx2ko_composition:
    """Figure S9C: transferred cell-type shares per genotype beside the shares published in Figure 5D of Zhu 2023."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx2ko_composition.py"),
        labels=data(KO + "transferred_labels.tsv.gz"),
        bars=data(KO + "zhu2023_figure5d_bar_segments.json"),
    output:
        by_label=result("figure5", "cdx2_knockout/genotype_by_label.tsv"),
        comparison=result("figure5", "cdx2_knockout/figure5d_comparison.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" labels="{input.labels}" bars="{input.bars}" '
        'by_label="{output.by_label}" comparison="{output.comparison}"'


rule figure5_cdx2ko_regulon_targets:
    """Figure S9D: knockout minus wild type on the targets of every reproducible regulon of Figure 5B in SC-SIC, each
    against its own 500 random gene sets; also the CDX2(+) regulons' own target sets (about 2 min)."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx2ko_regulon_targets.py"),
        symbol_updates=script("figure5", "target_symbol_updates.tsv"),
        counts=data(KO + "recounted_counts.h5ad"),
        labels=data(KO + "transferred_labels.tsv.gz"),
        groups=GENE_GROUPS,
        regulons=data(REG + "network_regulons.tsv.gz"),
        all_programs=result("figure5", "regulons/all_programs.tsv"),
        program_values=result("figure5", "regulons/program_values.tsv"),
    output:
        effects=result("figure5", "cdx2_knockout/regulon_targets/set_effects.tsv"),
        sets=result("figure5", "cdx2_knockout/regulon_targets/sets.tsv"),
        context=result("figure5", "cdx2_knockout/regulon_targets/context_regulons.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" counts="{input.counts}" labels="{input.labels}" groups="{input.groups}" '
        'regulons="{input.regulons}" all_programs="{input.all_programs}" program_values="{input.program_values}" '
        'symbol_updates="{input.symbol_updates}" '
        'effects="{output.effects}" sets="{output.sets}" context="{output.context}"'
