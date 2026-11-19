# How the deposited inputs were made

The workflow in this repository starts from the deposited intermediate files (the Zenodo deposit). This folder describes the steps that made those files from public data. Most of these steps ran on a computing cluster. The workflow does not rerun them.

Every public source is listed in `SOURCES.tsv`, with its accession, link and citation. The `upstream` column of `inputs/<figure>.tsv` says how each input file was made, and the sections below give the shared steps.

Where a detail is unknown, this page says so ("Unknown").

## (a) Reads and Cell Ranger recounts

- **Input:** raw reads (FASTQ or submitter BAM) from ENA and SRA for E-MTAB-8901, E-MTAB-9536, E-MTAB-9543, E-MTAB-10187, E-MTAB-10268, GSE125970, GSE139535 (submitter BAM files, converted to FASTQ with 10x bamtofastq), GSE148073, GSE151117, GSE167880, GSE185036, GSE185224, GSE197064, GSE199636 (multiome RNA), GSE202493, GSE202497, GSE214852, GSE221156 and GSE226346.
- **Tool:** Cell Ranger 9.0.1 with the 10x reference `refdata-gex-GRCh38-2024-A`, run per capture. Each recount was checked against the run lists in the public records.
  - The GSE151117 and GSE167880 grafts, the GSE185036 graft and the three GSE199636 grafts (human stem-cell-derived islets grafted into mice) were counted on the combined human and mouse reference `refdata-gex-GRCh38_and_GRCm39-2024-A`. Only their human genes are used.
  - The GSE214852 grafts (human intestinal organoids grown in mice) were counted on the human reference alone.
  - GSE214852 has three captures, each sequenced in two runs. They were counted run by run, and each capture was then rebuilt from its runs' `molecule_info.h5` files, counting each molecule (barcode, gene and UMI) once. These files come from the recounts in this step, not from the public deposits.
  - GSE185224 was counted with its gene expression and cell-hashing libraries together. Only the gene expression features are used.
  - GSE197064 (Olaniru 2023, fetal pancreas) was counted for all four of its pooled libraries; cells were assigned to specimens with the sample maps the authors deposited.
- **Output:** per-capture count matrices with 38,606 genes (the human genes of the reference).
- **Feeds:** `shared/ec_identity_pseudobulk/*`, `shared/sc_sic_sample_pseudobulk/*`, `shared/gse199636_sic_pseudobulk/*`, `shared/sc_endocrine_counts/*`, `shared/fitting_matrix/*`, `shared/gse214852_hio/*`, `shared/gut_reference_analysis.rds` (every study in it except Hickey 2023 and Egozi 2021), `figure4/pseudobulk/*` (except GSE114412), `figure4/gse199636_recount/*`, `figure4/graft_marker_counts/*`, `figure5/cdx2_knockout/recounted_counts.h5ad` and `figure_s6/olaniru2023_*`.
- **Not recounted:** these were used as deposited:
  - GSE114412: inDrops count tables;
  - Hickey 2023 (Dryad) and Egozi 2021 (Zenodo): count objects;
  - OMIX236 (Yu 2021b, fetal pancreas): the authors' count matrix.
- **Options:** `cellranger count` with 8 cores and 120 GB memory unless noted.
  - E-MTAB-9536 and E-MTAB-9543: `--include-introns=true --create-bam=false --chemistry=auto`. 35 of their 115 runs (10 of 26 and 25 of 89) used `--include-introns=true --create-bam=true` with the chemistry detected automatically (no `--chemistry` option).
  - E-MTAB-8901, E-MTAB-10187, E-MTAB-10268, GSE125970, GSE151117, GSE167880, GSE202497 and GSE214852: `--include-introns=true --create-bam=true`, chemistry detected automatically. The GSE151117 and GSE167880 grafts on the combined reference: `--include-introns=true --create-bam=false --chemistry=auto`.
  - GSE148073 and GSE221156: `--chemistry=auto --include-introns=true --create-bam=false --disable-ui`.
  - GSE139535: `--chemistry=SC3Pv2 --create-bam=false --disable-ui`, introns included (the Cell Ranger default), 16 cores and 184 GB.
  - GSE185036, GSE185224 and GSE226346: `--include-introns=true --create-bam=false --chemistry=auto --disable-ui`. GSE185036 also used `--r1-length=28`, and GSE185224 was counted with `--libraries` and `--feature-ref`.
  - GSE199636: `--chemistry=ARC-v1`, introns included, `--create-bam=false --disable-ui`.
  - GSE202493 (Zhu 2023 CDX2 wild type SRR19140209 and knockout SRR19140208, one library each): `--include-introns=true --create-bam=false`, chemistry detected automatically (Single Cell 3' v3), 16 cores and 120 GB.
  - **Unknown:** the `cellranger count` options of the GSE197064 runs.

## (b) Quality control and doublet calls

- **Input:** the count matrices from (a).
- **Tool:** one project-wide QC policy, and scDblFinder 1.24.10 under R 4.5.2. These ran in one pinned software environment.
  - Changing the thread count did not change the doublet calls.
- **Output:** QC-passed objects with a doublet call per cell.
- **Feeds:** the cell sets behind every file in (a). No figure uses the doublet call. The candidate doublets of Figure S1D are a separate flag based on marker read shares (`figure3_supp/s3/cells.tsv.gz` and `figure3_supp/s3/sc_sic_cells.tsv.gz`).
- **Thresholds:** the same rule for every recount, applied capture by capture. A cell passes when it has at least 200 genes detected and a mitochondrial share (genes whose symbol starts with MT-) of at most min(50%, max(20%, median + 3 × 1.4826 × MAD)) within its capture. scDblFinder ran per capture with `clusters = TRUE, dbr = NULL, dbr.sd = 1`, serially, with fixed seeds.
- **Use:** where a study's released cells are its authors' final post-QC set (the Gut Cell Atlas data and Hickey 2023), those cells are the analysis set. The figures use cells with published labels (c). The doublet call is kept as a flag and removes no cell. Hickey 2023 has no doublet call, because its complete captures cannot be identified in the deposited object.
- **GSE202493** (no published labels): the QC rule above defines its nuclei. Its scDblFinder call is kept as a column and removes no nucleus. 18,473 of 20,433 nuclei pass (wild type 7,478, knockout 10,995).

## (c) Attaching published labels

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
  - **GSE197064 (Olaniru 2023):** the study deposited no cell-level label table. Cell types were matched to the recounted cells from the authors' web viewer. 22,810 of the 28,879 cells received a cell type and a specimen; the others are left out.
  - **OMIX236 (Yu 2021b):** the cell annotation deposited with the authors' count matrix.
  - **Other studies:** the per-cell metadata in each study's GEO, Dryad or Zenodo deposit.
- **Harmonization:** study labels were mapped to a shared set:
  - intestinal: EC, other terminal enteroendocrine, unresolved and progenitor;
  - SC-islet: SC-SIC, beta, alpha, delta, epsilon and gamma. In Augsornworawat 2023 (GSE199636) the study's EB cells count as SC-SIC together with its EC cells.
  - Cells without a published label are left out of every figure.
- **Excluded donors:** the seven pediatric donors of Elmentaite 2020 (E-MTAB-8901) with Crohn's disease are left out of every analysis, including the fit of the four-marker EC classifier. Their disease status comes from the study's sample and data relationship file. `shared/gut_reference_roster/gut_reference_roster.json` lists them, and the workflow scripts read it.
- **Feeds:** every file whose description in `inputs/` says "published label", including `figure2/heldout_studies/*.h5ad` and `figure4/zhu2023_primary_labels.tsv.gz`.

## (d) SoupX ambient RNA correction

- **Input:** the raw (unfiltered) droplet matrices from (a) and the labels from (c).
- **Tool:** SoupX 1.6.2, run per 10x capture.
  - Soup profile: droplets with 1 to 100 UMIs.
  - Clusters: the published labels.
  - Contamination rate: `autoEstCont`, then `adjustCounts`.
- **Output:** raw and corrected counts of the cells kept.
- **Feeds:**
  - `figure3_supp/s3/gut_soupx_captures/*.rds`, one file per capture;
  - `figure3_supp/s3/sic_counts_raw_and_soupx.rds`;
  - `figure1/*_soupx_serotonin_genes.tsv.gz`;
  - `figure4/soupx_graft_*`;
  - `figure4/graft_marker_counts/soupx.tsv.gz` (section g).
- **Adult islet captures** (Fasolino 2022, Bandesh 2026): the clusters were the labeled endocrine cells only. One
  Fasolino 2022 capture (GSM4453624, 208 cells) is absent from `figure1/adult_islet_soupx_serotonin_genes.tsv.gz`
  because SoupX could not estimate its contamination.
- **Figure 4 graft captures:** SoupX could not estimate the contamination of one Balboa 2022 3-month graft capture
  (SRR13805905), so it enters `figure4/soupx_graft_*` and `figure4/graft_marker_counts/soupx.tsv.gz` uncorrected,
  within the pooled 3-month graft.

## (e) Pseudobulk sums

- **Input:** raw counts of labeled cells from (a) to (c).
- **Tool:** R and Python scripts that add up raw counts per study, sample or donor, and label. Each sum also records the number of cells with at least one count, and the cell and library totals.
- **Gene axes:**
  - intestinal: genes shared by the intestinal studies (21,520 genes);
  - SC-islet: the 38,563 genes of the 10x reference whose symbol occurs only once (38,606 genes minus 43 with a duplicated symbol).
  - The intestinal axis depends on the whole set of intestinal studies.
  - The Figure 4 trees (`figure4/pseudobulk/*`) keep each study's own genes: 38,606 for the recounts in (a), 16,540 for the GSE114412 inDrops tables.
- **Feeds:**
  - `shared/ec_identity_pseudobulk/*_sums.tsv.gz`, `*_detected.tsv.gz` and `*_groups.tsv`, and the fold changes and signatures computed from them (`shared/ec_identity/*`);
  - `shared/sc_sic_sample_pseudobulk/*` and `shared/gse199636_sic_pseudobulk/*`;
  - `figure4/pseudobulk/*` and `figure4/sc_alpha_profiles/*`.
- **Cell-level intestinal reference:** `shared/gut_reference_analysis.rds` holds the raw counts of the labeled intestinal cells on the same 21,520 genes, with their study, donor and label. Figures 4B, 4C, S1D and S8A read it.
- **Brain reference (Figure 2B):** made the same way from the CELLxGENE h5ad file of Siletti 2023 (clusters 397, serotonin neurons, and 395, dopamine neurons; single nuclei), summed per donor. The output is `figure3_supp/s9/siletti2023_neurons_*`. Figure 2B uses the serotonin neurons.
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
- **Scoring:** AUCell (`auc_threshold` 0.05, unnormalized, seed 906056) on each target dataset: the six SC-islet studies, the intestinal studies Elmentaite 2020 + 2021, Hickey 2023 and Egozi 2021, the two adult islet studies, and the McCauley 2023 HIO grafts.
  - Each regulon has 100 null regulons of the same size. Their targets are matched on mean count and detection fraction in 5 by 5 bins.
  - p value: empirical and two-sided, (1 + nulls at least as extreme) / 101.
  - A regulon is scored only if at least 10 of its targets, and at least half of them, are measured.
- **Output:** `figure5/regulons/network_regulons.tsv.gz`, `network_training_studies.tsv`, `regulon_target_effects.tsv.gz`, `mccauley2023_cells.tsv` and `mccauley2023_regulon_auc.tsv.gz`.

## (g) Inputs of the four-marker EC classifier

The label-free four-marker EC classifier is fitted inside the workflow (Figures S2 and S3). Its fitted constants are a result of the workflow, not a deposited file. These deposited files feed it and the panels that use it:

- **Fitting matrix** (`shared/fitting_matrix/*`): log1p counts per million (raw count divided by the cell's Gene Expression library total, times 10^6) of 1,954 enteroendocrine cells by 38,606 features. 1,199 are primary-tissue cells of the fitting studies (Yu 2021a, Elmentaite 2020, Wang 2020 and Burclaff 2022) and 755 are hPSC-derived HIO cells (E-MTAB-10268 and GSE214852). The cells come from the recounts in (a) and the labels in (c). The workflow leaves out the 35 cells of the excluded donors (c), which leaves 1,164 fitting cells.
- **Test studies** (`figure2/heldout_studies/*.h5ad`): the enteroendocrine cells of Elmentaite 2021 (897 cells), Hickey 2023 (3,533 nuclei) and Egozi 2021 (784 cells), with raw counts and their published labels grouped as EC, other terminal enteroendocrine, unresolved and progenitor. The HIO grafts of McCauley 2023 come from `shared/gse214852_hio/*`.
- **SC-islet cells:** `shared/sc_endocrine_counts/metadata.tsv`, `figure3_supp/s2/sc_panel_gene_counts.tsv.gz` (a gene subset of the SC endocrine count matrix), `figure3_supp/s2/sc_graphs.rds` (the per-study neighbor graphs, built with a lower-memory version of the classifier's graph builder that gives the same graphs) and `figure3/f3_03_scores_mature_endocrine.tsv.gz` (read for the cell order and published labels only).
- **Graft marker counts (Figure 4B):** `figure4/graft_marker_counts/` holds, for every SC-SIC and SC-beta cell (published labels) in the samples of the Figure 4B rows, the counts of the four markers TPH1, SLC18A1, FEV and LMX1A and the cell's library size.
  - `soupx.tsv.gz`: the eight graft rows, SoupX-corrected (d), with the column sum of the corrected matrix as library size (GSE151117, GSE167880, GSE139535, GSE185036 and GSE199636; 48,706 cells).
  - `uncorrected.tsv.gz`: the three extended-culture rows (GSE167880 stage 7 weeks 0 and 6, and the GSE199636 long-culture samples), raw counts with the total UMI count as library size (26,024 cells).
  - The workflow scores these cells with the fitted constants.

## (h) CDX2 knockout label transfer (GSE202493)

- **Input:**
  - reference: the 25,675 cells of Zhu 2023 (GSE202497) with published labels, from (a) to (c);
  - query: the GSE202493 nuclei recounted in (a) that pass QC in (b), 18,473 nuclei. Reference and query genes are matched by their Ensembl IDs (both GRCh38-2024-A, 38,606 genes).
- **Tool:** kNN label transfer with scanpy 1.10.0.
  - Settings: 2,000 reference highly variable genes, 50 principal components, 30 neighbors, confidence at least 0.6, seed 20260930.
- **Output:** `figure5/cdx2_knockout/transferred_labels.tsv.gz`, with `figure5/cdx2_knockout/recounted_counts.h5ad` (the recounted counts of the same nuclei).
  - The column `in_geo_deposit` marks the 12,615 nuclei that are also in the study's GEO deposit (`GSE202493_raw.h5ad.gz`); every one of them passes QC. Figures 5C, S9C and S9D use these 12,615 nuclei.
  - 15,399 of the 18,473 nuclei receive a label; the other 3,074 have too few genes to transfer. Five-fold cross-validated recall of SC-EC in the reference is 0.933.
- **Related:** `figure5/cdx2_knockout/zhu2023_figure5d_bar_segments.json` holds bar positions read from page 10 of the published article with PyMuPDF (Figure S9C compares them with the transferred labels).

## (i) Chromatin data

- **GSE199636 (Augsornworawat 2023):**
  - Peaks: the study's own per-sample cellranger-arc peak files for 10 samples, merged into 381,599 common hg38 peaks. Overlapping or book-ended peaks were merged, keeping chr1-22, X and Y. No new peak calling.
  - Counts: the study's fragment files (as submitted to GEO) were counted over those peaks for each nucleus, then summed per sample and published label. The genome sequence was 10x `refdata-gex-GRCh38-2024-A`.
  - Output: `figure5/chromatin/gse199636/*` (Figures 4D, 4E and S8C).
  - **Unknown:** the cellranger-arc version the authors used.
- **GSE202498 (Zhu 2023):**
  - Peaks: the study's own binary peak matrix, 65,255 nuclei by 462,765 hg19 peaks (500 bp), with nuclei named by the study's cluster names.
  - Counts: summed over 2 kb windows on either side of UCSC hg19 ncbiRefSeq transcription start sites. No new peak calling.
  - Output: `figure5/zhu_chromatin/promoter_counts_by_gene.tsv` and `groups.tsv` (Figure 4D).
- **Hickey 2023:**
  - Data: the study's snATAC fragments and its "Common Epithelial Peak Set" (509,414 hg38 peaks), from Dryad `doi:10.5061/dryad.0zpc8672f`.
  - Cells: 1,961 paired enteroendocrine nuclei, counted over those peaks and summed per sample and published label.
  - Output: `figure5/chromatin/hickey/*` (Figure S8B).
- **Gene coordinates:** `figure5/coordinates/genes.gtf.gz` is `genes/genes.gtf.gz` of the 10x reference package, unchanged.
- **Motif-family screen:**
  - Motifs: JASPAR 2024 CORE vertebrates non-redundant (879 motifs), matched in the peaks of the three datasets above with motifmatchr 1.32.0 and TFBSTools 1.48.0 under R 4.5.2. Zhu 2023 was run with the same scripts; its software versions are assumed to be the same.
  - Zhu 2023 motifs were matched on 9,974 selected peaks. Motif families group motifs with Jaccard similarity ≥ 0.4. The test is an exact conditional Mantel-Haenszel test, run in R 4.6.1.
  - Output: `figure5/chromatin_motif_screen_g2g3/*`, read for one annotation column of the Figure S9B tests.

## (j) MSigDB

- **Input:** `msigdb.2026.1.zip` from Zenodo record 18968178 (MD5 `512ba99c6827141a9d471972b812d4ac`, 44,335,790 bytes). This is the URL and checksum built into msigdbr 26.1.1.
- **Tool:** msigdbr 26.1.1, human gene sets (no cross-species mapping).
- **Use:** Hallmark and C8 (cell type) collections (Figures 2E and S5).
- **Deposit:** none. `workflow/scripts/common/fetch_url_inputs.py` downloads the archive and checks its SHA-256 against `inputs/` (`shared/msigdb/msigdb.2026.1.zip`). `msigdb_c8/` describes the C8 sets.
- **License:** MSigDB gene sets from v2022.1 on are CC BY 4.0. For any reported C8 set, also cite the source study named for that set.

## (k) Fetal pancreas pseudobulk (Figure S6)

- **Input:**
  - Olaniru 2023 (GSE197064): the recounts in (a) with the cell types in (c); 28,879 cells, 22,810 of them with a cell type and a specimen. Gene symbols that occur more than once are summed into one row.
  - Yu 2021b (OMIX236): the authors' count matrix and cell annotation of the whole-pancreas cells of seven embryos (30,422 cells, 22,550 genes); the endocrine-enriched sorts are left out.
- **Tool:** an R script that adds up the raw counts of every gene over the cells of each specimen (Olaniru 2023) or embryo (Yu 2021b) and published cell type. Every specimen and cell type is kept, including those with fewer than 10 cells, so the profile tables also give the label tallies. Gene rows keep the gene order of the source objects and profile columns keep a fixed order, because the matched random gene sets of Figure S6 depend on both orders.
- **Output:**
  - `figure_s6/olaniru2023_pseudobulk_counts.tsv.gz` (38,584 genes by 120 profiles) and `figure_s6/olaniru2023_profiles.tsv` (specimen, cell type, cells, age in weeks post-conception, 12 to 20);
  - `figure_s6/yu2021b_pseudobulk_counts.tsv.gz` (22,550 genes by 134 profiles) and `figure_s6/yu2021b_profiles.tsv` (embryo, cell type, cells, age in weeks as the authors give it, 8 to 19).
- The scoring of Figure S6 (profiles with at least 10 cells; Yu 2021b endocrine progenitors EP1 to EP4 pooled per embryo first) runs inside the workflow.

## (l) Display embedding (Figure 1C)

- **Input:** raw counts of the labeled endocrine cells of each study, from (a) to (c).
- **Tool:** per study, Seurat vst selection of 2,000 variable genes, log-normalization, scaling, 20 principal components, Harmony across samples (theta 2, lambda 1, sigma 0.1), then UMAP (uwot, cosine distance, 30 neighbors, min_dist 0.3, 500 epochs, seed 42).
- **Output:** `figure1/endocrine_display_cells.tsv.gz`: study, context, sample, published label, mapped class, UMAP coordinates and log-normalized TPH1, DDC, SLC18A1, CHGA, FEV, LMX1A, INS, GCG and SST for 139,701 cells. The embedding is for display only; no result is computed from it.

## (m) Published gene lists of two preprints (Figure 3C)

- **Input:** two supplementary files, used unchanged:
  - Liu 2025 (bioRxiv, doi 10.64898/2025.12.23.696311, version 1): Supplementary Table 5, gene programs (`Liu2025_SupplementaryTable5_gene_programs.csv`);
  - Fiancette 2026 (bioRxiv, doi 10.64898/2026.05.22.726825, version 1): Table S3, cluster marker genes of the H1 and H3 lines (`Fiancette2026_TableS3_cluster_markers.xlsx`).
- **Deposit:** none. The files are not redistributed; `workflow/scripts/common/fetch_url_inputs.py` downloads them from bioRxiv into `figure3/preprint_gene_lists/` and checks their SHA-256 against `inputs/figure3.tsv`. If the server refuses every attempt, the script says so; download the two files by hand from the URLs in `inputs/figure3.tsv` and put them at the paths given there.
