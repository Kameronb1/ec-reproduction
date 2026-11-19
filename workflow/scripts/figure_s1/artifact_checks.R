# Supplementary Figure S1D: do Figure 2C's SC-SIC vs intestinal EC differences survive ambient-RNA correction (SoupX) and removal of
# candidate doublets?
# Figure 2C's fit is rebuilt first and must reproduce the published Figure 2C table: one pseudobulk per intestinal EC donor
# (Elmentaite 2021, Hickey 2023, Egozi 2021, Elmentaite 2020 without its seven Crohn's disease donors; cells with published
# EC labels) and per SC-SIC sequenced sample (six SC studies), limma-voom on TMM-normalized counts of the
# filterByExpr genes (mitochondrial and ribosomal genes out), study as a duplicateCorrelation block with the Elmentaite pair
# one block (two voom / duplicateCorrelation rounds). It is then refitted, with the same tested genes, units and blocks, on
#   soupx     each unit's counts minus what SoupX removes from its cells (raw minus corrected counts per cell; intestinal
#             cells corrected per capture on every gene the capture measured; negative results clipped at zero);
#   doublets  each unit's counts without its cells with other-lineage reads (intestinal EC cells with > 2% of reads from enterocyte,
#             colonocyte or goblet markers; SC-SIC cells with > 5% from INS, IAPP, GCG, SST, PPY, GHRL); a unit left
#             without cells is dropped.
# Inputs (name=path):
#   analysis        analysis.rds: list(counts = genes x cells sparse raw counts, genes, metadata) of the primary intestinal
#                   cells (metadata: project_cell_id, family_id, biological_unit_id, author_group)
#   sc_groups       groups.tsv.gz of the per-sample SC pseudobulk (group_id, dataset, sample, population, cells,
#                   usable_feature_umi, file); sc_pb = folder holding its SC_SIC files (gene, count)
#   gse199636_groups groups.tsv of the GSE199636 SC-SIC per-sample pseudobulk; gse199636_pb = folder holding its files
#   donor_fit       contrast_donor_block_levels.tsv.gz (published Figure 2C table: gene, logFC, FDR, call)
#   gut_soupx       folder of the per-capture SoupX results (<capture>.rds: status, ec_raw, ec_soupx genes x EC cells)
#   sc_soupx        sic_counts_raw_and_soupx.rds: list(raw, soupx = genes x SC-SIC cells, cells = data.frame(cell, sample))
#   gut_doublets    cells.tsv.gz (cell, candidate_doublet) for intestinal EC cells
#   sc_doublets     sc_sic_cells.tsv.gz (cell, candidate_doublet) for SC-SIC cells
# Output:
#   out             genes_before_after.tsv.gz (per tested gene: Figure 2C call, logFC, FDR, and the same after each check)
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(Matrix); library(limma); library(edgeR) })
ok <- function(name, cond) if (!isTRUE(cond)) stop("check failed: ", name)

CROHNS_UNITS <- c("T017", "T019", "T176", "T189", "T197", "T202", "T203")   # Elmentaite 2020 child donors with Crohn's disease
MAIN <- c("GSE151117", "GSE167880", "GSE185036", "GSE199636", "GSE202497", "GSE226346")
FAMILIES <- c(GCA2021 = "SF_ELMENTAITE2021_GUT_ATLAS", HICKEY2023 = "SF_HICKEY2023_ADULT_EEC",
              EGOZI2021 = "SF_EGOZI2021_DEVELOPMENTAL_EEC", ELM2020 = "SF_ELMENTAITE2020_DEVGUT")
EXPECTED_EC <- c(GCA2021 = 319L, HICKEY2023 = 780L, EGOZI2021 = 170L, ELM2020 = 45L)
TECHNICAL <- "^MT-|^RP[SL][0-9]"
EC_GROUPS <- c("EXPLICIT_EC", "EXPLICIT_AUTHOR_EC")

# ---- Figure 2C's units -------------------------------------------------------------------------------------------------
primary <- readRDS(A[["analysis"]])
md <- primary$metadata
crohn <- md$family_id == FAMILIES[["ELM2020"]] & md$biological_unit_id %in% CROHNS_UNITS
ok("roster_35_cells", sum(crohn) == 35L)
primary$counts <- primary$counts[, !crohn, drop = FALSE]; md <- md[!crohn, , drop = FALSE]
rownames(primary$counts) <- primary$genes
for (k in names(FAMILIES)) ok(paste0("ec_cells_", k), sum(md$family_id == FAMILIES[[k]] & md$author_group %in% EC_GROUPS) == EXPECTED_EC[[k]])
nat_sym <- lapply(setNames(names(FAMILIES), names(FAMILIES)), function(k) {
  ec <- md$family_id == FAMILIES[[k]] & md$author_group %in% EC_GROUPS
  setNames(Matrix::rowSums(primary$counts[, ec, drop = FALSE]), primary$genes) })

groups <- unique(read.delim(gzfile(A[["sc_groups"]]), check.names = FALSE))
sic <- groups[groups$population == "SC_SIC" & groups$dataset %in% MAIN, ]
cnt <- lapply(seq_len(nrow(sic)), function(i) {
  a <- read.delim(gzfile(file.path(A[["sc_pb"]], sic$file[i]))); if (sum(a$count) != sic$usable_feature_umi[i]) stop("pseudobulk total")
  setNames(a$count, a$gene) })
ng <- read.delim(A[["gse199636_groups"]])
ncnt <- lapply(seq_len(nrow(ng)), function(i) { d <- read.delim(gzfile(file.path(A[["gse199636_pb"]], ng$file[i]))); setNames(d$count, d$gene) })
for (n in setdiff(names(sic), names(ng))) ng[[n]] <- NA
sic <- rbind(sic, ng[, names(sic)]); cnt <- c(cnt, ncnt)
pool_sc <- function(rows, genes) Reduce(`+`, lapply(rows, function(i) cnt[[i]][genes]))
genes_sc <- Reduce(intersect, lapply(which(sic$dataset %in% MAIN), function(i) names(cnt[[i]])))
COMMON <- sort(Reduce(intersect, c(list(genes_sc), lapply(nat_sym, names))))
COMMON <- COMMON[!grepl(TECHNICAL, COMMON, ignore.case = TRUE)]
nat_raw <- vapply(names(FAMILIES), function(k) as.numeric(nat_sym[[k]][COMMON]), numeric(length(COMMON)))
sc_raw <- vapply(MAIN, function(s) pool_sc(which(sic$dataset == s), COMMON), numeric(length(COMMON)))
rownames(nat_raw) <- rownames(sc_raw) <- COMMON
arm <- factor(c(rep("PRIMARY_EC", ncol(nat_raw)), rep("SC_SIC", ncol(sc_raw))), levels = c("PRIMARY_EC", "SC_SIC"))
dge <- edgeR::DGEList(counts = cbind(nat_raw, sc_raw), group = arm)
TREND_GENES <- rownames(dge)[edgeR::filterByExpr(dge, group = arm)]

# units: intestinal donors (cells listed so the corrections can act per cell), SC samples
UNITS <- list(); UCELLS <- list()
for (k in names(FAMILIES)) {
  ec <- md$family_id == FAMILIES[[k]] & md$author_group %in% EC_GROUPS
  for (u in sort(unique(md$biological_unit_id[ec]))) {
    UNITS[[length(UNITS) + 1L]] <- data.frame(arm = "PRIMARY_EC", study = k, unit = u)
    UCELLS[[length(UCELLS) + 1L]] <- md$project_cell_id[ec & md$biological_unit_id == u] } }
for (d in MAIN) for (smp in sort(unique(sic$sample[sic$dataset == d]))) {
  UNITS[[length(UNITS) + 1L]] <- data.frame(arm = "SC_SIC", study = d, unit = smp); UCELLS[[length(UCELLS) + 1L]] <- smp }
UNITS <- do.call(rbind, UNITS)
gut_prof <- function(cells) setNames(as.numeric(Matrix::rowSums(primary$counts[COMMON, cells, drop = FALSE])), COMMON)
prof0 <- lapply(seq_len(nrow(UNITS)), function(i) if (UNITS$arm[i] == "PRIMARY_EC") gut_prof(UCELLS[[i]]) else
  pool_sc(which(sic$dataset == UNITS$study[i] & sic$sample == UNITS$unit[i]), COMMON))
ok("units_44", nrow(UNITS) == 44L)

fit_units <- function(prof, keep_unit = rep(TRUE, nrow(UNITS))) {
  Yb <- do.call(cbind, prof[keep_unit]); colnames(Yb) <- paste(UNITS$study[keep_unit], UNITS$unit[keep_unit], sep = "|")
  armb <- factor(ifelse(UNITS$arm[keep_unit] == "SC_SIC", "SC_SIC", "PRIMARY_EC"), levels = c("PRIMARY_EC", "SC_SIC"))
  blk <- ifelse(UNITS$study[keep_unit] %in% c("GCA2021", "ELM2020"), "ELM_PAIR", UNITS$study[keep_unit])
  d <- edgeR::DGEList(counts = Yb, group = armb)
  d <- edgeR::calcNormFactors(d[rownames(d) %in% TREND_GENES, , keep.lib.sizes = FALSE], method = "TMM")
  des <- model.matrix(~ armb)
  v1 <- limma::voom(d, des); c1 <- limma::duplicateCorrelation(v1, des, block = blk)$consensus.correlation
  v2 <- limma::voom(d, des, block = blk, correlation = c1); c2 <- limma::duplicateCorrelation(v2, des, block = blk)$consensus.correlation
  f <- limma::eBayes(limma::lmFit(v2, des, block = blk, correlation = c2))
  tt <- limma::topTable(f, coef = "armbSC_SIC", number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, FDR = tt$adj.P.Val, call = ifelse(tt$adj.P.Val < 0.05 & tt$logFC > 0, "SC_SIC_HIGHER",
             ifelse(tt$adj.P.Val < 0.05 & tt$logFC < 0, "NATIVE_EC_HIGHER", "NO_CALL")), stringsAsFactors = FALSE)
}
D0 <- fit_units(prof0)
pub <- read.delim(gzfile(A[["donor_fit"]]), stringsAsFactors = FALSE)
D0 <- D0[match(pub$gene, D0$gene), ]
ok("refit_reproduces_fig2c", setequal(pub$gene, D0$gene) && nrow(pub) == nrow(D0) && identical(pub$call, D0$call) &&
   max(abs(pub$logFC - D0$logFC)) <= 1e-9 && max(abs(pub$FDR - D0$FDR) / pmax(pub$FDR, 1e-300)) <= 1e-6)

# ---- SoupX: subtract what SoupX removes from each unit's cells ---------------------------------------------------------
removed_on_common <- function(raw, soupx, cols) {
  g <- intersect(COMMON, rownames(raw)); R <- raw[g, cols, drop = FALSE] - soupx[g, cols, drop = FALSE]
  v <- setNames(numeric(length(COMMON)), COMMON); v[g] <- as.numeric(Matrix::rowSums(R)); v }
gut_ids <- unlist(UCELLS[UNITS$arm == "PRIMARY_EC"])
ok("gut_ec_cells_1314", length(gut_ids) == 1314L)
CK <- list.files(A[["gut_soupx"]], pattern = "[.]rds$", full.names = TRUE)
ok("captures_118", length(CK) == 118L)
REM <- Matrix::Matrix(0, length(COMMON), length(gut_ids), sparse = TRUE, dimnames = list(COMMON, gut_ids))
MEAS_N <- setNames(integer(length(gut_ids)), gut_ids)
for (f in CK) { x <- readRDS(f)
  ok(paste0("capture_corrected_", basename(f)), identical(x$status, "CORRECTED"))
  g <- intersect(COMMON, rownames(x$ec_raw)); cc <- colnames(x$ec_raw)
  REM[g, cc] <- x$ec_raw[g, , drop = FALSE] - x$ec_soupx[g, , drop = FALSE]
  MEAS_N[cc] <- length(g) }
ok("every_gut_cell_corrected", all(MEAS_N > 0) && setequal(names(MEAS_N), gut_ids))
ss <- readRDS(A[["sc_soupx"]])
prof_soupx <- prof0
for (i in seq_len(nrow(UNITS))) {
  cols <- if (UNITS$arm[i] == "PRIMARY_EC") UCELLS[[i]] else ss$cells$cell[ss$cells$sample == UNITS$unit[i]]
  if (UNITS$arm[i] == "SC_SIC") {
    ok(paste0("sc_cells_", UNITS$unit[i]), length(cols) == sum(sic$cells[sic$sample == UNITS$unit[i]]))
    g <- intersect(COMMON, rownames(ss$raw))
    ok(paste0("sc_raw_equals_pseudobulk_", UNITS$unit[i]), all(abs(Matrix::rowSums(ss$raw[g, cols, drop = FALSE]) - prof0[[i]][g]) == 0)) }
  rem <- if (UNITS$arm[i] == "PRIMARY_EC") setNames(as.numeric(Matrix::rowSums(REM[, cols, drop = FALSE])), COMMON) else removed_on_common(ss$raw, ss$soupx, cols)
  prof_soupx[[i]] <- pmax(prof0[[i]] - rem, 0) }
DS <- fit_units(prof_soupx)

# ---- doublets: drop the flagged cells ------------------------------------------------------------------------------------
fl_gut <- read.delim(gzfile(A[["gut_doublets"]]), stringsAsFactors = FALSE)
ok("doublet_flags_cover_gut_cells", all(gut_ids %in% fl_gut$cell))
fl_sc <- read.delim(gzfile(A[["sc_doublets"]]), stringsAsFactors = FALSE)
ok("doublet_flags_cover_sc_cells", all(ss$cells$cell %in% fl_sc$cell) && all(fl_sc$cell %in% ss$cells$cell))
BAD_GUT <- fl_gut$cell[fl_gut$candidate_doublet & fl_gut$cell %in% gut_ids]; BAD_SC <- fl_sc$cell[fl_sc$candidate_doublet]
prof_db <- prof0; keep_unit <- rep(TRUE, nrow(UNITS))
for (i in seq_len(nrow(UNITS))) {
  if (UNITS$arm[i] == "PRIMARY_EC") { keep <- setdiff(UCELLS[[i]], BAD_GUT)
    if (!length(keep)) { keep_unit[i] <- FALSE; next }
    prof_db[[i]] <- gut_prof(keep) } else {
    bad <- intersect(ss$cells$cell[ss$cells$sample == UNITS$unit[i]], BAD_SC)
    g <- intersect(COMMON, rownames(ss$raw)); v <- prof0[[i]]
    v[g] <- v[g] - as.numeric(Matrix::rowSums(ss$raw[g, bad, drop = FALSE])); prof_db[[i]] <- v } }
DD <- fit_units(prof_db, keep_unit)

# ---- per gene ------------------------------------------------------------------------------------------------------------
G <- data.frame(gene = pub$gene, call = pub$call, logFC = pub$logFC, FDR = pub$FDR,
                logFC_soupx = DS$logFC[match(pub$gene, DS$gene)], FDR_soupx = DS$FDR[match(pub$gene, DS$gene)], call_soupx = DS$call[match(pub$gene, DS$gene)],
                logFC_doublets = DD$logFC[match(pub$gene, DD$gene)], FDR_doublets = DD$FDR[match(pub$gene, DD$gene)], call_doublets = DD$call[match(pub$gene, DD$gene)])
G$core <- G$call != "NO_CALL" & abs(G$logFC) >= 1
ok("no_gene_lost_in_refits", !anyNA(G[, c("logFC_soupx", "logFC_doublets")]))
h <- gzfile(A[["out"]], "wt"); write.table(G, h, sep = "\t", quote = FALSE, row.names = FALSE); close(h)
