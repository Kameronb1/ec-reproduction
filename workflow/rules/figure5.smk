# Figure 5 (main, panels A-G) and its supplement (pages S17-S19): the tables behind each panel.
# Every rule takes the input check as its first input so nothing runs on unverified data.

F5_OK = str(RESULTS / "figure5_inputs_verified.txt")
PB = "shared/ec_identity_pseudobulk/"


rule figure5_regulator_rna:
    """Figure 5A (and, through rule figure5_cdx2_rna_contrast, the RNA axis of Figure S17E): thinned pseudobulk RNA of 34 regulators per cell type."""
    input:
        ok=F5_OK,
        script=script("figure5", "regulator_rna.R"),
        sc_sums=data(PB + "sc_sums.tsv.gz"), sc_groups=data(PB + "sc_groups.tsv"),
        gut_sums=data(PB + "gut_primary_sums.tsv.gz"), gut_groups=data(PB + "gut_primary_groups.tsv"),
        islet_sums=data(PB + "islet_sums.tsv.gz"), islet_groups=data(PB + "islet_groups.tsv"),
        induced_sums=data(PB + "induced_sums.tsv.gz"), induced_groups=data(PB + "induced_groups.tsv"),
        hio_sums=data(PB + "hio_mccauley_endocrine1_sums.tsv.gz"), hio_groups=data(PB + "hio_mccauley_endocrine1_groups.tsv"),
    output:
        by_sample=result("figure5", "regulator_rna/tf_expression_by_sample.tsv"),
        pooled=result("figure5", "regulator_rna/tf_expression_pooled_sparse_gut.tsv"),
        panel=result("figure5", "regulator_rna/f5_5a_regulator_rna.tsv"),
        hio_rows=result("figure5", "regulator_rna/f5_5a_regulator_rna__hio_mccauley_rows.tsv"),
        hio_detail=result("figure5", "regulator_rna/f5_s4_hio_mccauley_rna_detail.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" sc_sums="{input.sc_sums}" sc_groups="{input.sc_groups}" '
        'gut_primary_sums="{input.gut_sums}" gut_primary_groups="{input.gut_groups}" '
        'islet_sums="{input.islet_sums}" islet_groups="{input.islet_groups}" '
        'induced_sums="{input.induced_sums}" induced_groups="{input.induced_groups}" '
        'hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" '
        'by_sample="{output.by_sample}" pooled="{output.pooled}" panel="{output.panel}" '
        'hio_rows="{output.hio_rows}" hio_detail="{output.hio_detail}"'


rule figure5_cdx2ko_composition:
    """Figure S17F: transferred cell-type shares per genotype against Zhu 2023's published Figure 5D bars."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx2ko_composition.py"),
        labels=data("figure5/cdx2_knockout/transferred_labels.tsv.gz"),
        bars=data("figure5/cdx2_knockout/zhu2023_figure5d_bar_segments.json"),
    output:
        by_label=result("figure5", "cdx2_knockout/genotype_by_label.tsv"),
        comparison=result("figure5", "cdx2_knockout/figure5d_comparison.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" labels="{input.labels}" bars="{input.bars}" '
        'by_label="{output.by_label}" comparison="{output.comparison}"'


rule figure5_cdx2ko_group_test:
    """Figure 5F and S17G: knockout effect on gene groups in SC-SIC and SC-beta (and SC-alpha, scored after them
    and drawn in 5F) against 500 expression-matched random gene sets, and each set gene's own effect (per_gene.tsv,
    not drawn), on the Figure 3 gene groups G1-G3. Figure S17G draws Zhu 2023's Figure 5E genes from the same
    run (figure5e_genes.tsv)."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx2ko_group_test.py"),
        counts=data("figure5/cdx2_knockout/recounted_counts.h5ad"),
        labels=data("figure5/cdx2_knockout/transferred_labels.tsv.gz"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
    output:
        effects=result("figure5", "cdx2_knockout/group_test/group_effects.tsv"),
        subtracted=result("figure5", "cdx2_knockout/group_test/beta_subtracted.tsv"),
        genes=result("figure5", "cdx2_knockout/group_test/figure5e_genes.tsv"),
        nuclei=result("figure5", "cdx2_knockout/group_test/nuclei.tsv"),
        per_gene=result("figure5", "cdx2_knockout/group_test/per_gene.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" counts="{input.counts}" labels="{input.labels}" groups="{input.groups}" '
        'effects="{output.effects}" subtracted="{output.subtracted}" genes="{output.genes}" nuclei="{output.nuclei}" '
        'per_gene="{output.per_gene}"'


rule figure5_chromatin:
    """Figure 5C, 5D and S17A: promoter (and distal) accessibility of the gene groups in SC-islet snATAC. S17A uses the Figure 1
    labels for all three Augsornworawat 2023 samples (week 2 from its own Figure 1-label pseudobulks); 5D compares the
    grafts with the three week-2 samples Figure 4 uses, within SC-SIC, SC-beta and SC-alpha."""
    input:
        ok=F5_OK,
        script=script("figure5", "chromatin_promoters_culture_grafts.R"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        gtf=data("figure5/coordinates/genes.gtf.gz"),
        mtx=data("figure5/chromatin/gse199636/peaks_by_sample_label.mtx"),
        peaks=data("figure5/chromatin/gse199636/peak_ids.txt"),
        pbgroups=data("figure5/chromatin/gse199636/sample_label_groups.tsv"),
        w2_mtx=data("figure5/chromatin/gse199636_week2_fig1/peaks_by_sample_label.mtx"),
        w2_pbgroups=data("figure5/chromatin/gse199636_week2_fig1/sample_label_groups.tsv"),
        zhu_counts=data("figure5/zhu_chromatin/promoter_counts_by_gene.tsv"),
        zhu_groups=data("figure5/zhu_chromatin/groups.tsv"),
    output:
        promoters=result("figure5", "c5_promoters_by_sibling.tsv"),
        culture=result("figure5", "c5_culture_within_pops_wk2to4.tsv"),
        culture_genes=result("figure5", "c5_culture_per_gene.tsv"),
        culture_ext=result("figure5", "c5_culture_extended.tsv"),
        culture_ext_genes=result("figure5", "c5_culture_extended_per_gene.tsv"),
        graft=result("figure5", "c5_graft_relative.tsv"),
        graft_pop=result("figure5", "c5_graft_relative_by_population.tsv"),
        agreement=result("figure5", "c5_graft_window_agreement.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" groups="{input.groups}" gtf="{input.gtf}" mtx="{input.mtx}" peaks="{input.peaks}" '
        'pbgroups="{input.pbgroups}" w2_mtx="{input.w2_mtx}" w2_pbgroups="{input.w2_pbgroups}" zhu_counts="{input.zhu_counts}" zhu_groups="{input.zhu_groups}" '
        'out_promoters="{output.promoters}" out_culture="{output.culture}" out_culture_genes="{output.culture_genes}" '
        'out_culture_ext="{output.culture_ext}" out_culture_ext_genes="{output.culture_ext_genes}" '
        'out_graft="{output.graft}" out_graft_pop="{output.graft_pop}" out_agreement="{output.agreement}"'


rule figure5_chromatin_relative:
    """Figure 5B, S17B and S17C: 5C's measure (each gene's promoter accessibility against the typical gene of the same cell
    type, per gene group) in Zhu 2023 SC-islets on days 21, 32 and 39, Hickey 2023 native intestinal enteroendocrine cells
    in the small intestine and colon, and Singh 2024 induced enteroendocrine cells on days 5 and 6."""
    input:
        ok=F5_OK,
        script=script("figure5", "chromatin_promoters_relative.R"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        zhu_counts=data("figure5/zhu_chromatin/promoter_counts_by_gene.tsv"),
        zhu_groups=data("figure5/zhu_chromatin/groups.tsv"),
        gtf=data("figure5/coordinates/genes.gtf.gz"),
        hickey_mtx=data("figure5/chromatin/hickey/peaks_by_sample_label.mtx"),
        hickey_peaks=data("figure5/chromatin/hickey/peak_ids.txt"),
        hickey_pbgroups=data("figure5/chromatin/hickey/sample_label_groups.tsv"),
        refseq=data("figure5/coordinates/ncbiRefSeq_hg19.txt.gz"),
        singh_mtx=data("figure5/chromatin/gse238274/peaks_by_sample_label.mtx"),
        singh_peaks=data("figure5/chromatin/gse238274/peak_ids.txt"),
        singh_pbgroups=data("figure5/chromatin/gse238274/sample_label_groups.tsv"),
    output:
        relative=result("figure5", "c5_relative_zhu_hickey_singh.tsv"),
        per_gene=result("figure5", "c5_relative_zhu_hickey_singh_per_gene.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" groups="{input.groups}" zhu_counts="{input.zhu_counts}" zhu_groups="{input.zhu_groups}" '
        'gtf="{input.gtf}" hickey_mtx="{input.hickey_mtx}" hickey_peaks="{input.hickey_peaks}" hickey_pbgroups="{input.hickey_pbgroups}" '
        'refseq="{input.refseq}" singh_mtx="{input.singh_mtx}" singh_peaks="{input.singh_peaks}" singh_pbgroups="{input.singh_pbgroups}" '
        'out_groups="{output.relative}" out_genes="{output.per_gene}"'


rule figure5_identity_selection:
    input:
        ok=F5_OK,
        script=script("figure5", "identity_selection.py"),
        regulons=data("figure5/identity_factors/network_regulons.tsv.gz"),
        training=data("figure5/identity_factors/network_training_studies.tsv"),
        effects=data("figure5/identity_factors/regulon_target_effects.tsv.gz"),
        gene_bins=data("figure5/regulator_rna_targets_main/matched_null_gene_bins.tsv"),
        gene_groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        symbol_updates=script("figure5", "target_symbol_updates.tsv"),
    output:
        all_programs=result("figure5", "identity/all_programs.tsv"),
        factors=result("figure5", "identity/factors.tsv"),
        program_values=result("figure5", "identity/program_values.tsv"),
        gut_parts=result("figure5", "identity/gut_parts.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" regulons="{input.regulons}" training="{input.training}" effects="{input.effects}" '
        'gene_bins="{input.gene_bins}" gene_groups="{input.gene_groups}" symbol_updates="{input.symbol_updates}" '
        'all_programs="{output.all_programs}" '
        'factors="{output.factors}" program_values="{output.program_values}" gut_parts="{output.gut_parts}"'


rule figure5_identity_target_placement:
    input:
        ok=F5_OK,
        script=script("figure5", "identity_target_placement.py"),
        factors=result("figure5", "identity/factors.tsv"),
        all_programs=result("figure5", "identity/all_programs.tsv"),
        regulons=data("figure5/identity_factors/network_regulons.tsv.gz"),
        gene_bins=data("figure5/regulator_rna_targets_main/matched_null_gene_bins.tsv"),
        gene_groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        sums=data("shared/ec_identity_pseudobulk/sc_sums.tsv.gz"),
        groups=data("shared/ec_identity_pseudobulk/sc_groups.tsv"),
        motif_families=data("figure5/chromatin_motif_screen_g2g3/motif_family_definitions.tsv"),
        motif_hits=data("figure5/chromatin_motif_screen_g2g3/replicated_family_hits.tsv"),
        symbol_updates=script("figure5", "target_symbol_updates.tsv"),
    output:
        tests=result("figure5", "identity/tests.tsv"),
        factors_display=result("figure5", "identity/factors_display.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" factors="{input.factors}" all_programs="{input.all_programs}" regulons="{input.regulons}" '
        'gene_bins="{input.gene_bins}" gene_groups="{input.gene_groups}" pseudobulk_sums="{input.sums}" '
        'pseudobulk_groups="{input.groups}" motif_families="{input.motif_families}" motif_hits="{input.motif_hits}" '
        'symbol_updates="{input.symbol_updates}" tests="{output.tests}" factors_display="{output.factors_display}"'


rule figure5_identity_mccauley_check:
    input:
        ok=F5_OK,
        script=script("figure5", "identity_mccauley_check.py"),
        program_values=result("figure5", "identity/program_values.tsv"),
        effects=data("figure5/identity_factors/regulon_target_effects.tsv.gz"),
        training=data("figure5/identity_factors/network_training_studies.tsv"),
        cells=data("figure5/identity_factors/mccauley2023_cells.tsv"),
        auc=data("figure5/identity_factors/mccauley2023_regulon_auc.tsv.gz"),
    output:
        out=result("figure5", "identity/gse214852_values_primary_labels.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" program_values="{input.program_values}" effects="{input.effects}" '
        'training="{input.training}" cells="{input.cells}" auc="{input.auc}" out="{output.out}"'


rule figure5_identity_supplement_values:
    input:
        ok=F5_OK,
        script=script("figure5", "identity_supplement_values.R"),
        values=result("figure5", "identity/program_values.tsv"),
        parts=result("figure5", "identity/gut_parts.tsv"),
        factors=result("figure5", "identity/factors.tsv"),
        tests=result("figure5", "identity/tests.tsv"),
        fdisp=result("figure5", "identity/factors_display.tsv"),
    output:
        out=result("figure5", "identity/f5_s_gh_values_ALL_ELIGIBLE.tsv"),
        singh=result("figure5", "identity/f5_s_gh_singh_values_ALL_ELIGIBLE.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" values="{input.values}" parts="{input.parts}" factors="{input.factors}" '
        'tests="{input.tests}" fdisp="{input.fdisp}" out="{output.out}" singh_out="{output.singh}"'


rule figure5_cdx_sites:
    """Figure S17E, chromatin axis: CDX motif sites by study (the table also keeps each EC population's CDX2 RNA). CDX2 RNA comes from the regulator RNA table made by
    rule figure5_regulator_rna (tf_expression_by_sample.tsv); McCauley 2023's RNA-only row is that rule's
    f5_5a_regulator_rna__hio_mccauley_rows.tsv and is not recomputed here."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx_sites_vs_rna.R"),
        family=data("figure5/chromatin_mechanism/motif_family_vs_rna.tsv"),
        families=data("figure5/chromatin_mechanism/motif_families.tsv"),
        zhu=data("figure5/zhu_chromatin/zhu_candidate_tf_directions.tsv"),
        rna=result("figure5", "regulator_rna/tf_expression_by_sample.tsv"),
    output:
        out=result("figure5", "cdx_sites/f5_s4_cdx_sites_vs_rna.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" family="{input.family}" families="{input.families}" zhu="{input.zhu}" '
        'rna="{input.rna}" out="{output.out}"'


rule figure5_cdx2_rna_contrast:
    """Figure S17E, RNA axis: CDX2 RNA, EC minus the study's other endocrine cells, in the same five studies as the
    chromatin axis. Profiles that Figure 5A draws keep rule figure5_regulator_rna's values (tf_expression_by_sample.tsv);
    the others are drawn here with the same depth, draws and seed."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx2_rna_ec_minus_other.R"),
        sc_sums=data(PB + "sc_sums.tsv.gz"), sc_groups=data(PB + "sc_groups.tsv"),
        gut_sums=data(PB + "gut_primary_sums.tsv.gz"), gut_groups=data(PB + "gut_primary_groups.tsv"),
        induced_sums=data(PB + "induced_sums.tsv.gz"), induced_groups=data(PB + "induced_groups.tsv"),
        hio_sums=data(PB + "hio_mccauley_endocrine1_sums.tsv.gz"), hio_groups=data(PB + "hio_mccauley_endocrine1_groups.tsv"),
        by_sample=result("figure5", "regulator_rna/tf_expression_by_sample.tsv"),
    output:
        out=result("figure5", "cdx_sites/f5_s4_cdx2_rna_ec_minus_other.tsv"),
        by_profile=result("figure5", "cdx_sites/f5_s4_cdx2_rna_by_profile.tsv"),
    shell:
        '"{RSCRIPT}" "{input.script}" sc_sums="{input.sc_sums}" sc_groups="{input.sc_groups}" '
        'gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" '
        'induced_sums="{input.induced_sums}" induced_groups="{input.induced_groups}" '
        'hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" by_sample="{input.by_sample}" '
        'out="{output.out}" by_profile="{output.by_profile}"'


rule figure5_cdx2ko_regulon_targets:
    """Figure S17H and S17I: the CDX2 knockout's effect on the targets of CDX2's own regulons (H; SC-SIC, SC-beta and
    SC-alpha, pseudobulk and per-cell detection) and on every scorable Figure 5E regulon in SC-SIC (I), each against its
    own 500 expression-matched random gene sets, with the readouts and null of rule figure5_cdx2ko_group_test. Targets named
    by an older symbol are read under their current one (target_symbol_updates.tsv; FAM49B is CYRIB)."""
    input:
        ok=F5_OK,
        script=script("figure5", "cdx2ko_regulon_targets.py"),
        symbol_updates=script("figure5", "target_symbol_updates.tsv"),
        counts=data("figure5/cdx2_knockout/recounted_counts.h5ad"),
        labels=data("figure5/cdx2_knockout/transferred_labels.tsv.gz"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        regulons=data("figure5/identity_factors/network_regulons.tsv.gz"),
        all_programs=result("figure5", "identity/all_programs.tsv"),
        program_values=result("figure5", "identity/program_values.tsv"),
    output:
        effects=result("figure5", "cdx2_knockout/regulon_targets/set_effects.tsv"),
        subtracted=result("figure5", "cdx2_knockout/regulon_targets/beta_subtracted.tsv"),
        sets=result("figure5", "cdx2_knockout/regulon_targets/sets.tsv"),
        context=result("figure5", "cdx2_knockout/regulon_targets/context_regulons.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" counts="{input.counts}" labels="{input.labels}" groups="{input.groups}" '
        'regulons="{input.regulons}" all_programs="{input.all_programs}" program_values="{input.program_values}" '
        'symbol_updates="{input.symbol_updates}" '
        'effects="{output.effects}" subtracted="{output.subtracted}" sets="{output.sets}" context="{output.context}"'
