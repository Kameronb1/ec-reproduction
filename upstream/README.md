# How the deposited inputs were made

The workflow in this repository starts from the deposited intermediate files (the Zenodo deposit). This folder describes the steps that made those files from public data. These steps ran on a computing cluster. The workflow does not rerun them.

Every public source is listed in `SOURCES.tsv`, with its accession, link and citation. The `upstream` column of `inputs/<figure>.tsv` says which step below made each deposited file.

Where a detail is unknown, this page says so ("Unknown").

## (a) Reads and Cell Ranger recounts

- **Input:** raw reads (FASTQ or submitter BAM) from ENA and SRA for E-MTAB-8901, E-MTAB-9536, E-MTAB-9543, E-MTAB-10187, E-MTAB-10268, GSE125970, GSE139535 (submitter BAM files, converted to FASTQ with 10x bamtofastq), GSE148073, GSE151117, GSE167880, GSE185036, GSE185224, GSE199636 (multiome RNA), GSE202493, GSE202497, GSE214852, GSE221156, GSE226346 and GSE238274.
- **Tool:** Cell Ranger 9.0.1 with the 10x reference `refdata-gex-GRCh38-2024-A`, run per capture. Each recount was checked against the run lists in the public records.
  - The GSE151117 and GSE167880 grafts, the GSE185036 graft and the three GSE199636 grafts (human stem-cell-derived islets grafted into mice) were counted on the combined human and mouse reference `refdata-gex-GRCh38_and_GRCm39-2024-A`. Only their human genes are used.
  - The GSE214852 grafts (human intestinal organoids grown in mice) were counted on the human reference alone.
  - GSE214852 (three captures, each sequenced in two runs) and day 6 of GSE238274 (one capture sequenced in two runs) were counted run by run. Each capture was then rebuilt from its runs' `molecule_info.h5` files, counting each molecule (barcode, gene and UMI) once. These files come from the recounts in this step, not from the public deposits.
  - GSE185224 was counted with its gene expression and cell-hashing libraries together. Only the gene expression features are used.
- **Output:** per-capture count matrices with 38,606 genes (the human genes of the reference).
- **Feeds:** `shared/ec_identity_pseudobulk/*`, `shared/sc_sic_sample_pseudobulk/*`, `shared/gse199636_sic_pseudobulk/*`, `shared/sc_endocrine_counts/*`, `shared/fitting_matrix/*`, `shared/gse214852_hio/*`, `shared/gut_reference_analysis.rds` (every study in it except Hickey 2023 and Egozi 2021), `figure2/heldout_studies/gse238274.h5ad`, `figure2/classifier_inputs/secondary_validation.rds`, `figure3/h_test_sources_four_marker/*`, `figure4/pseudobulk/*`, `figure4/gse199636_recount/*`, `figure5/cdx2_natural_variation/sc_islet_exports/*`, `figure5/gse199636_endpoint_objects/*`, `figure5/cdx2_knockout/recounted_counts.h5ad`.
- **Not recounted:** these were used as deposited:
  - GSE114412: inDrops count tables;
  - Hickey 2023 (Dryad) and Egozi 2021 (Zenodo): count objects.
- **Options:** `cellranger count` with 8 cores and 120 GB memory unless noted.
  - E-MTAB-9536 and E-MTAB-9543: `--include-introns=true --create-bam=false --chemistry=auto`. 35 of their 115 runs (10 of 26 and 25 of 89) used `--include-introns=true --create-bam=true` with the chemistry detected automatically (no `--chemistry` option).
  - E-MTAB-8901, E-MTAB-10187, E-MTAB-10268, GSE125970, GSE151117, GSE167880, GSE202497 and GSE214852: `--include-introns=true --create-bam=true`, chemistry detected automatically. The GSE151117 and GSE167880 grafts on the combined reference: `--include-introns=true --create-bam=false --chemistry=auto`.
  - GSE238274: `--chemistry=SC3Pv3 --include-introns=true --create-bam=true`.
  - GSE148073 and GSE221156: `--chemistry=auto --include-introns=true --create-bam=false --disable-ui`.
  - GSE139535: `--chemistry=SC3Pv2 --create-bam=false --disable-ui`, introns included (the Cell Ranger default), 16 cores and 184 GB.
  - GSE185036, GSE185224 and GSE226346: `--include-introns=true --create-bam=false --chemistry=auto --disable-ui`. GSE185036 also used `--r1-length=28`, and GSE185224 was counted with `--libraries` and `--feature-ref`.
  - GSE199636: `--chemistry=ARC-v1`, introns included, `--create-bam=false --disable-ui`.
  - GSE202493 (Zhu 2023 CDX2 wild type SRR19140209 and knockout SRR19140208, one library each): `--include-introns=true --create-bam=false`, chemistry detected automatically (Single Cell 3' v3), 16 cores and 120 GB.

## (b) Quality control and doublet calls

- **Input:** the count matrices from (a).
- **Tool:** one project-wide QC policy, and scDblFinder 1.24.10 under R 4.5.2. These ran in one pinned software environment.
  - Changing the thread count did not change the doublet calls.
- **Output:** QC-passed objects with a doublet call per cell.
- **Feeds:** the cell sets behind every file in (a). No figure uses the doublet call. The candidate doublets of Figure S9 are a separate flag based on marker read shares (`figure3_supp/s3/cells.tsv.gz` and `figure3_supp/s3/sc_sic_cells.tsv.gz`).
- **Thresholds:** the same rule for every recount, applied capture by capture. A cell passes when it has at least 200 genes detected and a mitochondrial share (genes whose symbol starts with MT-) of at most min(50%, max(20%, median + 3 × 1.4826 × MAD)) within its capture. scDblFinder ran per capture with `clusters = TRUE, dbr = NULL, dbr.sd = 1`, serially, with fixed seeds.
- **Use:** where a study's released cells are its authors' final post-QC set (the Gut Cell Atlas data, Hickey 2023, Singh 2024), those cells are the analysis set. The figures use cells with primary study labels (c). The doublet call is kept as a flag and removes no cell. Hickey 2023 has no doublet call, because its complete captures cannot be identified in the deposited object.
- **GSE202493** (no primary study labels): the QC rule above defines its nuclei. Its scDblFinder call is kept as a column and removes no nucleus. 18,473 of 20,433 nuclei pass (wild type 7,478, knockout 10,995).

## (c) Attaching primary study labels

- **Input:** the per-cell annotations each study published, matched to the recounted cells by exact capture and barcode.
- **Sources:**
  - **E-MTAB-8901:** the two CELLxGENE h5ad files (fetal and pediatric).
    - EC and other enteroendocrine labels come from an exact cell-key crosswalk to the later Gut Cell Atlas h5ad. This is a later marker-informed reannotation of the same cells, not a new label.
    - "EC cells (TAC1+)" become EC. D, I, K, L, M/X and N cells become other enteroendocrine cells.
    - 19 cells had no match and are left out.
  - **E-MTAB-9536 and E-MTAB-9543:** the Gut Cell Atlas h5ad.
  - **E-MTAB-10187:** the authors' Mendeley Data files and CELLxGENE.
  - **GSE151117:** a label file supplied by the authors. It is not publicly deposited.
    - The join matched 8,309 cells exactly.
  - **GSE139535 (WS4):** a per-cell cluster file supplied by the first author. It is not publicly deposited.
  - **Other studies:** the per-cell metadata in each study's GEO, Dryad or Zenodo deposit.
- **Harmonization:** study labels were mapped to a shared set:
  - intestinal: EC, other terminal enteroendocrine, unresolved and progenitor;
  - SC-islet: SC-SIC, beta, alpha, delta, epsilon and gamma. In Augsornworawat 2023 (GSE199636) the study's EB cells count as SC-SIC together with its EC cells.
  - Cells without a primary study label are left out of every figure.
- **Feeds:** every file whose description says "primary study label", including `figure2/heldout_studies/*.h5ad`, `figure3/primary_label_population_census.tsv` and `figure4/zhu2023_primary_labels.tsv.gz`.

## (d) SoupX ambient RNA correction

- **Input:** the raw (unfiltered) droplet matrices from (a) and the labels from (c).
- **Tool:** SoupX 1.6.2, run per 10x capture.
  - Soup profile: droplets with 1 to 100 UMIs.
  - Clusters: the primary study labels.
  - Contamination rate: `autoEstCont`, then `adjustCounts`.
- **Output:** raw and corrected counts of the cells kept.
- **Feeds:**
  - `figure3_supp/s3/gut_soupx_captures/*.rds`, one file per capture;
  - `figure3_supp/s3/sic_counts_raw_and_soupx.rds`;
  - `figure1/*_soupx_serotonin_genes.tsv.gz`;
  - `figure4/soupx_graft_*`;
  - `figure4/ec_calls_per_cell_soupx.tsv.gz`.
- **Adult islet captures** (Fasolino 2022, Bandesh 2026): the clusters were the labeled endocrine cells only. One
  Fasolino 2022 capture (GSM4453624, 208 cells) is absent from `figure1/adult_islet_soupx_serotonin_genes.tsv.gz`
  because SoupX could not estimate its contamination.
- **Figure 4 graft captures:** SoupX could not estimate the contamination of one Balboa 2022 3-month graft capture
  (SRR13805905), so it enters `figure4/soupx_graft_*` and `figure4/ec_calls_per_cell_soupx.tsv.gz` uncorrected, within
  the pooled 3-month graft.

## (e) Pseudobulk sums

- **Input:** raw counts of labeled cells from (a) to (c).
- **Tool:** R and Python scripts that add up raw counts per study, sample or donor, and label. Each sum also records the number of cells with at least one count, and the cell and library totals.
- **Gene axes:**
  - intestinal: genes shared by the intestinal studies (21,520 genes);
  - SC-islet: the 38,563 genes of the 10x reference whose symbol occurs only once (38,606 genes minus 43 with a duplicated symbol).
  - The intestinal axis depends on the whole set of intestinal studies.
  - The Figure 4 trees (`figure4/pseudobulk/*`) keep each study's own genes: 38,606 for the recounts in (a), 16,540 for the GSE114412 inDrops tables.
- **Feeds:**
  - `shared/ec_identity_pseudobulk/*_sums.tsv.gz`, `*_detected.tsv.gz` and `*_groups.tsv`;
  - `shared/sc_sic_sample_pseudobulk/*` and `shared/gse199636_sic_pseudobulk/*`;
  - `figure4/pseudobulk/*` and `figure4/sc_alpha_profiles/*`;
  - `figure3/h_test_sources_four_marker/four_marker_abstain_*`: the intestinal test sources Yu 2021, Wang 2020 (GSE125970) and Burclaff 2022 (GSE185224), summed per donor and by the four-marker EC classifier's result for each cell (classified as EC, abstained, or not classified as EC). Figure 3H reads them, and so do Figure 4F's ticks. The 10% detection floor that Figure 3H and Figure 4F apply to the gene list is computed inside the workflow.
- **Brain reference (Figure S12):** made the same way from the CELLxGENE h5ad files of Siletti 2023 (clusters 397 and 395, single nuclei) and Braun 2023 (medulla cluster 52, whole cells). The output is `figure3_supp/s9/*`.
- **Unknown:** the exact commands of the intestinal pseudobulk step. Its output files are checked against their checksums in `inputs/<figure>.tsv`.

## (f) Regulon networks (SCENIC) and regulon scores

- **Input:** per-cell expression of labeled cells.
- **Networks:** one gene regulatory network per study, 12 in all. Each network is named after the study it was built from, except `support_small_adult_intestine`, which was built from E-MTAB-9543.
  - GSE151117, GSE167880, GSE185036, GSE199636, GSE202497 and GSE226346;
  - E-MTAB-9536 and E-MTAB-9543;
  - Hickey 2023 and Egozi 2021;
  - GSE148073 and GSE221156.
- **Tool:** pySCENIC 0.12.1 with arboreto 0.1.6 and ctxcore 0.2.0 (Python 3.10.6), in the `aertslab/pyscenic` container pinned by digest. Each regulon is a factor's predicted target set.
  - Network step: GRNBoost2.
  - Motif step: the hg38 v10_clust cisTarget ranking databases (10 kb up and down of the TSS; 500 bp up and 100 bp down), the `motifs-v10nr_clust` annotation and the `allTFs_hg38` factor list.
  - Runs: five per network (seeds 907051 to 907055). Each run uses 80% of cells, at most 250 per sample and population.
  - Consensus: a regulon is kept if it appears in at least 3 of the 5 runs and has at least 10 targets. Activating and repressing regulons are kept separate.
  - The GSE148073 network was trained without donor HPAP034.
  - The GSE199636 network used these settings. It was trained on SC-SIC (the study's EC and EB nuclei) and SC-beta nuclei of SC-islet 1 and 2, on the study's Figure 1 labels (12,312 nuclei). Each run drew 250 nuclei per sample and population (1,000 nuclei) and kept genes detected in at least 10 of them (1%), leaving out mitochondrial (MT-) and ribosomal protein (RPS, RPL) genes (17,364 to 17,485 genes per run). Commands: `pyscenic grn --method grnboost2 --seed <run seed>`, then `pyscenic ctx --all_modules` with both ranking databases and the motif annotation. 44 regulons passed the consensus rule.
- **Scoring:** AUCell (`auc_threshold` 0.05, unnormalized, seed 906056) on each target dataset.
  - Each regulon has 100 null regulons of the same size. Their targets are matched on mean count and detection fraction in 5 by 5 bins.
  - p value: empirical and two-sided, (1 + nulls at least as extreme) / 101.
  - A regulon is scored only if at least 10 of its targets, and at least half of them, are measured.
- **Output:**
  - `figure5/identity_factors/network_regulons.tsv.gz`, `network_training_studies.tsv`, `regulon_target_effects.tsv.gz` and `mccauley2023_regulon_auc.tsv.gz`;
  - `figure4/regulon_screen.tsv.gz` and `figure4/nfe2l1_regulon_unit_activity.tsv`.

## (g) Figure 2 classifier run and matched-random test

- **Input:** two training objects:
  - `figure2/classifier_inputs/secondary_validation.rds`;
  - `shared/gut_reference_analysis.rds`, built from the intestinal studies' count matrices and primary study labels.
- **Tool:** a leave-one-study-out run of the 20-gene classifier.
  - It ran in a pinned R 4.2.0 container with limma, Matrix 1.6.1.1 and single-threaded BLAS.
  - Seven gene-ranking datasets from six studies were used. Singh 2024 was scored as an additional test study, with the same thresholds.
  - The same run made the label-swap null replicates.
- **Output:** `figure2/classifier_run/*`, `shared/figure2_panel_d/D_effect_heatmap.tsv` and `figure2_supp/classifier_run/recurrence_null_replicates.tsv`.
- **Matched-random test:** 200 random gene sets per study.
  - It writes `figure2_supp/classifier_run/matched_random.tsv` and `hio_matched_random.tsv`.
  - Its script is in this repository (`workflow/scripts/figure2_supp/classifier.R`, section 2). It takes about an hour on one core, so the default workflow reads the deposited tables instead of rerunning it.
- **Unknown:** the limma version inside the container.

## (h) CDX2 knockout label transfer (GSE202493)

- **Input:**
  - reference: the 25,675 cells of Zhu 2023 (GSE202497) with primary study labels, from (a) to (c);
  - query: the GSE202493 nuclei recounted in (a) that pass QC in (b), 18,473 nuclei. Reference and query genes are matched by their Ensembl IDs (both GRCh38-2024-A, 38,606 genes).
- **Tool:** kNN label transfer with scanpy 1.10.0.
  - Settings: 2,000 reference highly variable genes, 50 principal components, 30 neighbors, confidence at least 0.6, seed 20260930.
- **Output:** `figure5/cdx2_knockout/transferred_labels.tsv.gz`, with `figure5/cdx2_knockout/recounted_counts.h5ad` (the recounted counts of the same nuclei).
  - The column `in_geo_deposit` marks the 12,615 nuclei that are also in the study's GEO deposit (`GSE202493_raw.h5ad.gz`); every one of them passes QC. Figures 5F and S17F to S17I use these 12,615 nuclei.
  - 15,399 of the 18,473 nuclei receive a label; the other 3,074 have too few genes to transfer. Five-fold cross-validated recall of SC-EC in the reference is 0.933.
- **Related:** `figure5/cdx2_knockout/zhu2023_figure5d_bar_segments.json` holds bar positions read from page 10 of the published article with PyMuPDF.

## (i) Chromatin data

- **GSE199636 (Augsornworawat 2023):**
  - Peaks: the study's own per-sample cellranger-arc peak files for 10 samples, merged into 381,599 common hg38 peaks. Overlapping or book-ended peaks were merged, keeping chr1-22, X and Y. No new peak calling.
  - Counts: the study's fragment files (as submitted to GEO) were counted over those peaks for each nucleus, then summed per sample and primary study label. The genome sequence was 10x `refdata-gex-GRCh38-2024-A`.
  - Output: `figure5/chromatin/gse199636/*`.
  - **Unknown:** the cellranger-arc version the authors used.
- **GSE202498 (Zhu 2023):**
  - Peaks: the study's own binary peak matrix, 65,255 nuclei by 462,765 hg19 peaks (500 bp), with nuclei named by the study's cluster names.
  - Counts: summed over 2 kb windows on either side of UCSC hg19 ncbiRefSeq transcription start sites. No new peak calling.
  - Output: `figure5/zhu_chromatin/promoter_counts_by_gene.tsv` and `groups.tsv`.
- **Hickey 2023:**
  - Data: the study's snATAC fragments and its "Common Epithelial Peak Set" (509,414 hg38 peaks), from Dryad `doi:10.5061/dryad.0zpc8672f`.
  - Cells: 1,961 paired enteroendocrine nuclei.
- **Singh 2024:**
  - Data: GEO GSE238275 (normal-line days 1, 3, 5 and 6).
  - Peaks were called by this project: MACS3 3.0.4 per primary study label (`-f BEDPE -g hs --nomodel --keep-dup all -q 0.01 --call-summits`), 501 bp around each summit, then merged into 186,413 hg19 peaks.
- **Motif activity (chromVAR):**
  - Motifs: JASPAR 2024 CORE vertebrates non-redundant (879 motifs).
  - Software: motifmatchr 1.32.0, TFBSTools 1.48.0, chromVAR 1.32.0 (50 background iterations, seed 20260925), R 4.5.2.
  - Applied to the four datasets above. Zhu 2023 was run with the same scripts; its software versions are assumed to be the same.
  - Outputs: `figure5/zhu_chromatin/zhu_candidate_tf_directions.tsv`, `figure5/chromatin_mechanism/*` and `figure5/chromatin_motif_screen_g2g3/*`.
- **Motif-family screen:** reuses those motif matches. Zhu 2023 motifs were matched on 9,974 selected peaks. Motif families group motifs with Jaccard similarity ≥ 0.4. The test is an exact conditional Mantel-Haenszel test, run in R 4.6.1.

## (j) MSigDB

- **Input:** `msigdb.2026.1.zip` from Zenodo record 18968178 (MD5 `512ba99c6827141a9d471972b812d4ac`, 44,335,790 bytes). This is the URL and checksum built into msigdbr 26.1.1.
- **Tool:** msigdbr 26.1.1, human gene sets (no cross-species mapping).
- **Use:** Hallmark and C8 (cell type) collections.
- **Deposit:** none. The workflow downloads the archive directly and checks its SHA-256 against `inputs/` (`shared/msigdb/msigdb.2026.1.zip`).
- **License:** MSigDB gene sets from v2022.1 on are CC BY 4.0. For any reported C8 set, also cite the source study named for that set.
