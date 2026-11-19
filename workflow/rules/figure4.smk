# Figure 4: SC-islet grafts and extended culture. A, mature endocrine composition of each graft against its in vitro
# samples; B, cells classified as EC by the four-marker classifier and the four markers' expression change; C, gap to
# primary intestinal EC on the gene groups; D, E, promoter and distal chromatin accessibility of the gene groups. The
# composition table of Supplementary Figure S7 and the chromatin tables of Supplementary Figure S8B-C are made here too.
# The four-marker constants come from the classifier unit (figure_s2_s3). Rules marked resources heavy_r=1 hold the
# intestinal reference matrix or a chromatin peak matrix in memory; with --resources heavy_r=1 they run one at a time.

FIG4_OK = str(RESULTS / "figure4_inputs_verified.txt")
F4_TREES = ["GSE167880", "GSE151117", "GSE139535", "GSE185036_S6D2", "GSE185036_10wpt", "GSE202497"]


def f4_tree_files(n):
    """Inputs tree1..tree<n> (profiles, genes, counts) in the order the scripts expect."""
    files = {}
    for i, t in enumerate(F4_TREES[:n], start=1):
        for part, name in [("profiles", "profiles.tsv"), ("genes", "genes.tsv"), ("counts", "counts.tsv.gz")]:
            files[f"tree{i}_{part}"] = data(f"figure4/pseudobulk/{t}/{name}")
    return files


def f4_args(wildcards, input):
    """name="path" for every named input except the input-check file (ok) and the script."""
    return " ".join(f'{k}="{v}"' for k, v in input.items() if k not in ("ok", "script"))


rule figure4_composition:
    """Figure 4A and Figure S7: mature endocrine composition of every SC-islet sample (published labels)."""
    input:
        ok=FIG4_OK,
        script=script("figure4", "composition.R"),
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
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{input.script}" {params.args} f4b="{output.f4b}" s1="{output.s1}"'


rule figure4_ec_calls:
    """Figure 4A (samples per row end) and 4B (four-marker scores and calls of every SC-SIC and SC-beta cell)."""
    input:
        ok=FIG4_OK,
        script=script("figure4", "ec_calls.R"),
        markers_soupx=data("figure4/graft_marker_counts/soupx.tsv.gz"),
        markers_uncorrected=data("figure4/graft_marker_counts/uncorrected.tsv.gz"),
        constants=result("figure_s2_s3", "fitted_constants.json"),
        comparisons=data("figure4/comparisons_whole_cell.tsv"),
        comparisons_nuclei=data("figure4/comparisons_nuclei.tsv"),
    output:
        rows=result("figure4", "row_samples.tsv"),
        cells=result("figure4", "C_cells.tsv.gz"),
        shares=result("figure4", "four_marker_row_shares.tsv"),
    params:
        args=f4_args,
    shell:
        '"{RSCRIPT}" "{input.script}" {params.args} rows="{output.rows}" cells="{output.cells}" shares="{output.shares}"'


rule figure4_levels:
    """Figure 4B (the four markers' expression change) and 4C (gap to primary intestinal EC on the gene groups)."""
    input:
        ok=FIG4_OK,
        script=script("figure4", "levels.R"),
        **f4_tree_files(6),
        gut=data("shared/gut_reference_analysis.rds"),
        cfap92=data("figure4/gut_cfap92_alias_row.rds"),
        universe=data("shared/gene_groups/universe_paneld_tested.txt"),
        nuc_profiles=data("figure4/gse199636_recount/profiles.tsv"),
        nuc_pseudobulk=data("figure4/gse199636_recount/pseudobulk_counts.tsv.gz"),
        soupx_profiles=data("figure4/soupx_graft_profiles.tsv"),
        soupx_pseudobulk=data("figure4/soupx_graft_pseudobulk.tsv.gz"),
        rows=result("figure4", "row_samples.tsv"),
        groups=data("shared/fig3_gene_groups.tsv.gz"),
        roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
    output:
        bars=result("figure4", "C_bars.tsv"),
        gaps=result("figure4", "E_group_gaps.tsv"),
        end_spread=result("figure4", "E_end_spread.tsv"),
    params:
        args=f4_args,
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{input.script}" {params.args} bars="{output.bars}" gaps="{output.gaps}" end_spread="{output.end_spread}"'


rule figure4_chromatin:
    """Figure 4D (culture weeks 2-4, months 6 and 12, graft median), 4E and S8C (grafts against random gene sets):
    promoter and distal accessibility of the gene groups in Augsornworawat 2023."""
    input:
        ok=FIG4_OK,
        script=script("figure4", "chromatin_promoters_culture_grafts.R"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        gtf=data("figure5/coordinates/genes.gtf.gz"),
        mtx=data("figure5/chromatin/gse199636/peaks_by_sample_label.mtx"),
        peaks=data("figure5/chromatin/gse199636/peak_ids.txt"),
        pbgroups=data("figure5/chromatin/gse199636/sample_label_groups.tsv"),
    output:
        culture=result("figure4", "c5_culture_within_pops_wk2to4.tsv"),
        culture_genes=result("figure4", "c5_culture_per_gene.tsv"),
        culture_ext=result("figure4", "c5_culture_extended.tsv"),
        culture_ext_genes=result("figure4", "c5_culture_extended_per_gene.tsv"),
        graft=result("figure4", "c5_graft_relative.tsv"),
        graft_pop=result("figure4", "c5_graft_relative_by_population.tsv"),
        agreement=result("figure4", "c5_graft_window_agreement.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{input.script}" groups="{input.groups}" gtf="{input.gtf}" mtx="{input.mtx}" peaks="{input.peaks}" '
        'pbgroups="{input.pbgroups}" out_culture="{output.culture}" out_culture_genes="{output.culture_genes}" '
        'out_culture_ext="{output.culture_ext}" out_culture_ext_genes="{output.culture_ext_genes}" '
        'out_graft="{output.graft}" out_graft_pop="{output.graft_pop}" out_agreement="{output.agreement}"'


rule figure4_chromatin_relative:
    """Figure 4D (Zhu 2023, days 21, 32 and 39) and S8B (Hickey 2023, small intestine and colon): each gene group's
    promoter accessibility against the typical gene of the same cell type."""
    input:
        ok=FIG4_OK,
        script=script("figure4", "chromatin_promoters_relative.R"),
        groups=data("shared/gene_groups/gene_groups_main/gene_groups_definitions.tsv.gz"),
        zhu_counts=data("figure5/zhu_chromatin/promoter_counts_by_gene.tsv"),
        zhu_groups=data("figure5/zhu_chromatin/groups.tsv"),
        gtf=data("figure5/coordinates/genes.gtf.gz"),
        hickey_mtx=data("figure5/chromatin/hickey/peaks_by_sample_label.mtx"),
        hickey_peaks=data("figure5/chromatin/hickey/peak_ids.txt"),
        hickey_pbgroups=data("figure5/chromatin/hickey/sample_label_groups.tsv"),
    output:
        relative=result("figure4", "c5_relative_zhu_hickey.tsv"),
        per_gene=result("figure4", "c5_relative_zhu_hickey_per_gene.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{input.script}" groups="{input.groups}" zhu_counts="{input.zhu_counts}" zhu_groups="{input.zhu_groups}" '
        'gtf="{input.gtf}" hickey_mtx="{input.hickey_mtx}" hickey_peaks="{input.hickey_peaks}" hickey_pbgroups="{input.hickey_pbgroups}" '
        'out_groups="{output.relative}" out_genes="{output.per_gene}"'
