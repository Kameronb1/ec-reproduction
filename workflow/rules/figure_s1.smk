# Supplementary Figure S1: controls for Figures 1 and 2. A: islet hormone transcripts per cell group (its table is
# written by the Figure 1 unit, hormones_by_group.tsv); B: Figure 2A's agreement with intestinal EC when every donor or
# sample counts once; C: each population's fold changes against the intestinal EC signature, gene by gene; D: Figure
# 2C's genes after ambient-RNA correction and after removing cells with other cell types' reads.
# Rules marked resources heavy_r=1 load the 2.2 GB SC-islet count matrix or the intestinal cell matrix; with
# --resources heavy_r=1 they run one at a time.

S1_OK = str(RESULTS / "figure_s1_inputs_verified.txt")
S1_PB = "shared/ec_identity_pseudobulk/"
S1_PB_FILES = [S1_PB + t + s for t in ("gut_primary", "sc", "hio_mccauley_endocrine1")
               for s in ("_sums.tsv.gz", "_detected.tsv.gz", "_groups.tsv")]
S1_SC_COUNTS = ["shared/sc_endocrine_counts/" + f for f in
                ("metadata.tsv", "genes.txt", "counts_indptr.i32", "counts_indices.i32", "counts_data.i32")]
_S1_INPUTS = [r["file"] for r in read_tsv(REPO / "inputs" / "figure_s1.tsv")]
S1_SC_SAMPLES = [f for f in _S1_INPUTS if f.startswith("shared/sc_sic_sample_pseudobulk/")]
S1_GSE199636 = [f for f in _S1_INPUTS if f.startswith("shared/gse199636_sic_pseudobulk/")]
S1_SOUPX_CAPTURES = [f for f in _S1_INPUTS if f.startswith("figure3_supp/s3/gut_soupx_captures/")]


rule figure_s1_agreement_per_donor:
    input:
        ok=S1_OK,
        script=script("figure_s1", "agreement_per_donor.py"),
        lib=script("figure_s1", "signature_lib.py"),
        pb=[data(f) for f in S1_PB_FILES],
        sc=[data(f) for f in S1_SC_COUNTS],
    output:
        result("figure_s1", "dots_pooled_vs_balanced.tsv"),
    params:
        pb=data("shared/ec_identity_pseudobulk"),
        sc=data("shared/sc_endocrine_counts"),
    resources:
        heavy_r=1,
    shell:
        '"{PYTHON}" "{input.script}" lib="{input.lib}" pb="{params.pb}" sc_counts="{params.sc}" out="{output}"'


rule figure_s1_signature_scatter:
    input:
        ok=S1_OK,
        script=script("figure_s1", "signature_scatter.py"),
        lib=script("figure_s1", "signature_lib.py"),
        pb=[data(S1_PB + "hio_mccauley_endocrine1" + s) for s in ("_sums.tsv.gz", "_detected.tsv.gz", "_groups.tsv")],
        gut=data("shared/ec_identity/gut_source_fold_changes.tsv.gz"),
        signatures=data("shared/ec_identity/q1_gene_signatures.tsv.gz"),
    output:
        genes=result("figure_s1", "gene_signatures.tsv.gz"),
        concordance=result("figure_s1", "concordance.tsv"),
    params:
        pb=data("shared/ec_identity_pseudobulk"),
    shell:
        '"{PYTHON}" "{input.script}" lib="{input.lib}" pb="{params.pb}" gut="{input.gut}" signatures="{input.signatures}" '
        'out_genes="{output.genes}" out_concordance="{output.concordance}"'


rule figure_s1_share_two_fold:
    input:
        ok=S1_OK,
        script=script("figure_s1", "share_two_fold.py"),
        genes=result("figure_s1", "gene_signatures.tsv.gz"),
    output:
        result("figure_s1", "share_two_fold.tsv"),
    shell:
        '"{PYTHON}" "{input.script}" genes="{input.genes}" out="{output}"'


rule figure_s1_artifact_checks:
    input:
        ok=S1_OK,
        script=script("figure_s1", "artifact_checks.R"),
        analysis=data("shared/gut_reference_analysis.rds"),
        sc_files=[data(f) for f in S1_SC_SAMPLES],
        g6_files=[data(f) for f in S1_GSE199636],
        donor_fit=data("shared/sc_sic_vs_gut_donor_block/contrast_donor_block_levels.tsv.gz"),
        captures=[data(f) for f in S1_SOUPX_CAPTURES],
        sc_soupx=data("figure3_supp/s3/sic_counts_raw_and_soupx.rds"),
        gut_doublets=data("figure3_supp/s3/cells.tsv.gz"),
        sc_doublets=data("figure3_supp/s3/sc_sic_cells.tsv.gz"),
    output:
        result("figure_s1", "genes_before_after.tsv.gz"),
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
        'gse199636_groups="{params.g6_groups}" gse199636_pb="{params.g6_pb}" donor_fit="{input.donor_fit}" '
        'gut_soupx="{params.captures}" sc_soupx="{input.sc_soupx}" gut_doublets="{input.gut_doublets}" '
        'sc_doublets="{input.sc_doublets}" out="{output}"'
