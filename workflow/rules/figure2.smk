# Figure 2: tables behind panels A-E. Scripts are in workflow/scripts/figure2/. The gene-group results (gene_groups.tsv.gz,
# mixed_genes.tsv, universe.txt) and the donor fit are read by other figures.
# R model fits carry resources heavy_r=1; run with --resources heavy_r=1 to keep them one at a time.

F2_OK = str(RESULTS / "figure2_inputs_verified.txt")
F2_S = script("figure2", "").rstrip("/\\")
F2_AGG = "shared/ec_identity_pseudobulk/"
F2_SC_SAMPLES = [
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


def f2_agg(key, stem):
    """The three files of one pseudobulk, as key_sums= key_detected= key_groups= arguments."""
    return {f"{key}_sums": data(F2_AGG + stem + "_sums.tsv.gz"),
            f"{key}_detected": data(F2_AGG + stem + "_detected.tsv.gz"),
            f"{key}_groups": data(F2_AGG + stem + "_groups.tsv")}


def f2_args(named):
    return " ".join(f'{k}="{v}"' for k, v in named.items())


F2_A_IN = {**f2_agg("gut", "gut_primary"), **f2_agg("sc", "sc"), **f2_agg("hio", "hio_mccauley_endocrine1"),
           "roster": data("shared/gut_reference_roster/gut_reference_roster.json")}

rule figure2_signature_correlations:
    input:
        ok=F2_OK,
        **F2_A_IN,
    output:
        dots=result("figure2", "signature_dots.tsv"),
        hio_dots=result("figure2", "signature_hio_dots.tsv"),
        pairwise=result("figure2", "signature_pairwise.tsv"),
        hio_pairwise=result("figure2", "signature_hio_pairwise.tsv"),
        gut_folds=result("figure2", "gut_ec_fold_changes.tsv.gz"),
    params:
        args=f2_args(F2_A_IN),
    shell:
        '"{PYTHON}" "{F2_S}/signature_correlations.py" {params.args} dots="{output.dots}" hio_dots="{output.hio_dots}" '
        'pairwise="{output.pairwise}" hio_pairwise="{output.hio_pairwise}" gut_folds="{output.gut_folds}"'


# Figure 2B: only the populations the panel draws enter the normalization and gene selection (dopamine neurons left out).
rule figure2_brain_similarity:
    input:
        ok=F2_OK,
        gut_sums=data(F2_AGG + "gut_primary_sums.tsv.gz"),
        gut_groups=data(F2_AGG + "gut_primary_groups.tsv"),
        sc_sums=data(F2_AGG + "sc_sums.tsv.gz"),
        sc_groups=data(F2_AGG + "sc_groups.tsv"),
        hio_sums=data(F2_AGG + "hio_mccauley_endocrine1_sums.tsv.gz"),
        hio_groups=data(F2_AGG + "hio_mccauley_endocrine1_groups.tsv"),
        islet_sums=data(F2_AGG + "islet_sums.tsv.gz"),
        islet_groups=data(F2_AGG + "islet_groups.tsv"),
        gut_roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
        islet_roster=data("shared/adult_islet_roster/roster.json"),
        brain_sums=data("figure3_supp/s9/siletti2023_neurons_sums.tsv.gz"),
        brain_groups=data("figure3_supp/s9/siletti2023_neurons_groups.tsv"),
    output:
        similarity=result("figure2", "brain_similarity.tsv"),
        unrounded=result("figure2", "brain_similarity_unrounded.tsv"),
    shell:
        '"{RSCRIPT}" "{F2_S}/brain_similarity.R" gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" '
        'sc_sums="{input.sc_sums}" sc_groups="{input.sc_groups}" hio_sums="{input.hio_sums}" hio_groups="{input.hio_groups}" '
        'islet_sums="{input.islet_sums}" islet_groups="{input.islet_groups}" gut_roster="{input.gut_roster}" '
        'islet_roster="{input.islet_roster}" brain_sums="{input.brain_sums}" brain_groups="{input.brain_groups}" '
        'drop="Brain dopamine neurons" similarity="{output.similarity}" similarity_unrounded="{output.unrounded}"'


rule figure2_donor_fit:
    input:
        ok=F2_OK,
        gut_sums=data(F2_AGG + "gut_primary_sums.tsv.gz"),
        gut_groups=data(F2_AGG + "gut_primary_groups.tsv"),
        roster=data("shared/gut_reference_roster/gut_reference_roster.json"),
        samples=[data("shared/sc_sic_sample_pseudobulk/" + f) for f in F2_SC_SAMPLES],
    output:
        levels=result("figure2", "donor_fit.tsv.gz"),
        profiles=result("figure2", "pooled_voom_profiles.tsv.gz"),
        tested=result("figure2", "tested_genes.txt"),
    params:
        sc_dir=data("shared/sc_sic_sample_pseudobulk"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{F2_S}/donor_fit.R" gut_sums="{input.gut_sums}" gut_groups="{input.gut_groups}" '
        'roster="{input.roster}" sc_dir="{params.sc_dir}" levels="{output.levels}" '
        'profiles="{output.profiles}" tested="{output.tested}"'


F2_D_IN = {**f2_agg("sc", "sc"), "gut_sums": data(F2_AGG + "gut_primary_sums.tsv.gz"),
           "gut_groups": data(F2_AGG + "gut_primary_groups.tsv"), "roster": data("shared/gut_reference_roster/gut_reference_roster.json")}

rule figure2_gene_groups:
    input:
        ok=F2_OK,
        gut_folds=result("figure2", "gut_ec_fold_changes.tsv.gz"),
        fit=result("figure2", "donor_fit.tsv.gz"),
        profiles=result("figure2", "pooled_voom_profiles.tsv.gz"),
        **F2_D_IN,
    output:
        definitions=result("figure2", "gene_groups.tsv.gz"),
        mixed=result("figure2", "mixed_genes.tsv"),
        universe=result("figure2", "universe.txt"),
    params:
        args=f2_args(F2_D_IN),
    shell:
        '"{PYTHON}" "{F2_S}/gene_groups.py" {params.args} gut_folds="{input.gut_folds}" fit="{input.fit}" '
        'profiles="{input.profiles}" definitions="{output.definitions}" mixed="{output.mixed}" universe="{output.universe}"'


rule figure2_enrich_pathways:
    input:
        ok=F2_OK,
        groups=result("figure2", "gene_groups.tsv.gz"),
        universe=result("figure2", "universe.txt"),
        fit=result("figure2", "donor_fit.tsv.gz"),
        msigdb=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all_terms=result("figure2", "pathways_all_terms.tsv.gz"),
        called=result("figure2", "pathways_called.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{F2_S}/enrich_pathways.R" groups="{input.groups}" universe="{input.universe}" '
        'fit="{input.fit}" msigdb_zip="{input.msigdb}" all_terms="{output.all_terms}" called="{output.called}"'


rule figure2_enrich_c8:
    input:
        ok=F2_OK,
        groups=result("figure2", "gene_groups.tsv.gz"),
        universe=result("figure2", "universe.txt"),
        fit=result("figure2", "donor_fit.tsv.gz"),
        msigdb=data("shared/msigdb/msigdb.2026.1.zip"),
    output:
        all_terms=result("figure2", "celltypes_all_terms.tsv"),
        called=result("figure2", "celltypes_called.tsv"),
    resources:
        heavy_r=1,
    shell:
        '"{RSCRIPT}" "{F2_S}/enrich_c8.R" groups="{input.groups}" universe="{input.universe}" '
        'fit="{input.fit}" msigdb_zip="{input.msigdb}" all_terms="{output.all_terms}" called="{output.called}"'
