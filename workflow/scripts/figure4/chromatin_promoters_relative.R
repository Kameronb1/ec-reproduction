# Figure 4D (Zhu 2023) and Supplementary Figure S8B (Hickey 2023): the promoter measure of
# chromatin_promoters_culture_grafts.R in two more single-nucleus ATAC studies. Within one cell type at one time point (or
# intestinal region), each gene's promoter accessibility (accessible events in peaks within 2 kb of any TSS of the gene,
# log2(1 + CPM)) minus the median of every gene with a promoter peak; per gene group (G1-G3), the median gene and the 25th
# and 75th percentiles of the group's genes. Each cell type is first binomially thinned to its shallowest time point or
# region: 20 seeded draws per cell type and point, each gene's log2(1 + CPM) averaged over the draws before the median
# and quartiles are taken.
#
#   4D    Zhu 2023 (GSE202498) SC-islets, days 21, 32 and 39 of differentiation (samples MM157, MM166, MM168): SC-SIC,
#         SC-beta and SC-alpha (the published clusters SCEC, SCbeta, SCalpha). Per-gene promoter counts (peaks within
#         2 kb of any UCSC hg19 NCBI RefSeq TSS of the gene).
#   S8B   Hickey 2023 primary adult intestine: EC and the other enteroendocrine cells (D, I, K, L, Mo and S cells,
#         pooled), small intestine (duodenum, jejunum, ileum) and colon (ascending to sigmoid), every donor pooled. hg38
#         peaks; promoters from GENCODE v44 transcript TSS, as in chromatin_promoters_culture_grafts.R.
#
# Inputs (name=path):
#   groups      gene group table (gene_groups_definitions.tsv.gz; columns group, gene)
#   zhu_counts  Zhu 2023 promoter counts per gene and group (promoter_counts_by_gene.tsv)
#   zhu_groups  Zhu 2023 group nuclei and total accessible events (groups.tsv)
#   gtf         GENCODE v44 / 10x GRCh38-2024-A genes.gtf.gz (transcript TSS, Hickey 2023)
#   hickey_mtx, hickey_peaks, hickey_pbgroups  Hickey 2023 peak x pseudobulk counts, peak ids (chr:start-end), pseudobulk columns
# Outputs: out_groups (one row per study, cell type, day or region, and gene group), out_genes (each group gene's value)
args <- commandArgs(trailingOnly = TRUE); a <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
suppressPackageStartupMessages({library(digest); library(Matrix); library(GenomicRanges); library(data.table)})
PROM <- 2000L; SEED <- 20261103L; REPS <- 20L
wt <- function(x, p) fwrite(x, p, sep = "\t", quote = TRUE, na = "NA")
lcpm <- function(x, depth) log2(1 + x / depth * 1e6)
seed_of <- function(key) SEED + strtoi(substr(digest(key, algo = "xxhash32", serialize = FALSE), 1, 7), 16L)
thin <- function(x, d) if (sum(x) > d) rbinom(length(x), x, d / sum(x)) else x
# mean over REPS draws of one seeded stream; level_of(draw) returns one draw's per-gene levels (same genes every draw)
avg_draws <- function(key, level_of) { set.seed(seed_of(key)); Reduce(`+`, lapply(seq_len(REPS), function(i) level_of())) / REPS }
GG <- as.data.table(read.delim(gzfile(a[["groups"]]), check.names = FALSE, stringsAsFactors = FALSE))
GG <- GG[, .(gene, g = c(G1_EC_GENES_SC_SIC_LACK = "G1", G2_EC_GENES_SC_SIC_HAVE = "G2", G3_SC_SIC_OWN_PROGRAM = "G3")[group])]
stopifnot(!anyNA(GG$g))
peaks_of <- function(ids) { p <- do.call(rbind, strsplit(sub(":", "-", ids), "-")); GRanges(p[, 1], IRanges(as.integer(p[, 2]), as.integer(p[, 3]))) }
promoter_map <- function(peaks, chrom, tss, gene) {
  po <- findOverlaps(peaks, resize(GRanges(chrom, IRanges(tss, width = 1)), 2 * PROM + 1L, fix = "center"))
  unique(data.table(gene = gene[subjectHits(po)], peak = queryHits(po)))
}
read_pb <- function(mtx, peaks, pbgroups) {
  G <- fread(pbgroups); M <- as(readMM(mtx), "CsparseMatrix"); ids <- readLines(peaks)
  stopifnot(nrow(M) == length(ids), ncol(M) == nrow(G))
  list(G = G, M = M, peaks = peaks_of(ids))
}

# one cell type at one point: lv holds every gene with a promoter peak (gene, level)
relative <- function(lv, key, depth, nuclei, donors = NA_integer_) {
  allm <- median(lv$level)
  genes <- merge(lv[, .(gene, level)], GG, by = "gene")[, `:=`(all_genes_median = allm, rel = level - allm)]
  groups <- genes[, .(genes = .N, rel = median(rel), q25 = quantile(rel, .25), q75 = quantile(rel, .75)), by = g][
    , `:=`(nuclei = nuclei, donors = donors, depth = depth, all_genes = nrow(lv))]
  list(groups = cbind(key, groups), genes = cbind(key, genes))
}
pt <- function(study, population, day = NA_integer_, region = NA_character_)
  data.table(study = study, population = population, day = as.integer(day), region = region)

# ---- 4D, Zhu 2023: per-gene promoter counts; thinning to the cell type's shallowest day, on total accessible events ----
zp <- fread(a[["zhu_counts"]]); zg <- fread(a[["zhu_groups"]])
ZDAY <- c(MM157 = 21L, MM166 = 32L, MM168 = 39L); ZPOP <- c(SCEC = "SC-SIC", SCbeta = "SC-beta", SCalpha = "SC-alpha")
RZ <- unlist(lapply(names(ZPOP), function(cl) {
  z <- zg[cluster == cl & sample %in% names(ZDAY)]; stopifnot(nrow(z) == length(ZDAY))
  tot <- setNames(z$total_accessible_events, z$sample); d <- min(tot)
  lapply(names(ZDAY), function(s) {
    x0 <- zp[[paste0(cl, "|", s)]]
    lv <- data.table(gene = zp$gene, level = avg_draws(paste("zhu", cl, s),
      function() lcpm(if (tot[[s]] > d) rbinom(length(x0), x0, d / tot[[s]]) else x0, d)))
    relative(lv, pt("Zhu 2023", ZPOP[[cl]], day = ZDAY[[s]]), d, z[sample == s, nuclei])
  })
}), recursive = FALSE)

# ---- S8B, Hickey 2023: EC and the other enteroendocrine cells, small intestine and colon ------------------------------
gl <- readLines(gzfile(a[["gtf"]])); f <- strsplit(gl[grepl("\ttranscript\t", gl, fixed = TRUE)], "\t", fixed = TRUE); rm(gl)
TX <- data.table(chrom = vapply(f, `[`, "", 1), start = as.integer(vapply(f, `[`, "", 4)), end = as.integer(vapply(f, `[`, "", 5)), strand = vapply(f, `[`, "", 7),
                 gene = sub('.*gene_name "([^"]+)".*', "\\1", vapply(f, `[`, "", 9))); rm(f)
TX[, tss := ifelse(strand == "+", start, end)]; TX <- unique(TX[, .(chrom, tss, gene)])
H <- read_pb(a[["hickey_mtx"]], a[["hickey_peaks"]], a[["hickey_pbgroups"]]); HMAP <- promoter_map(H$peaks, TX$chrom, TX$tss, TX$gene)
REG <- list(`small intestine` = c("Duodenum", "Proximal-jejunum", "Mid-jejunum", "Ileum"), colon = c("Ascending", "Transverse", "Descending", "Sigmoid"))
stopifnot(all(unique(H$G$region) %in% unlist(REG)))
HPOP <- list(EC = function(l) l == "Enterochromaffin", `other enteroendocrine` = function(l) l != "Enterochromaffin")
HKEY <- c(EC = "EC", `other enteroendocrine` = "other entero-\nendocrine", `small intestine` = "small\nintestine", colon = "colon")   # seed keys
col_h <- function(r, p) which(H$G$scope == "primary" & H$G$aggregation_level == "sample" & H$G$region %in% REG[[r]] & HPOP[[p]](H$G$analysis_label))
RH <- unlist(lapply(names(HPOP), function(p) {
  X <- lapply(setNames(names(REG), names(REG)), function(r) as.numeric(rowSums(H$M[, col_h(r, p), drop = FALSE]))); d <- min(vapply(X, sum, 0))
  lapply(names(REG), function(r) {
    lv <- HMAP[, .(n = sum(X[[r]][peak])), by = gene][, n := NULL]
    lv[, level := avg_draws(paste("hickey", HKEY[[p]], HKEY[[r]]), function() { x <- thin(X[[r]], d); HMAP[, .(n = sum(x[peak])), by = gene][lv, on = "gene", lcpm(n, d)] })]
    i <- col_h(r, p); relative(lv, pt("Hickey 2023", p, region = r), d, sum(H$G$cells[i]), uniqueN(H$G$donor[i]))
  })
}), recursive = FALSE)
rm(H)

# rows in the order computed (study, cell type, day or region), then by group and gene
R <- c(RZ, RH)
grp <- rbindlist(lapply(R, `[[`, "groups")); gen <- rbindlist(lapply(R, `[[`, "genes"))
o <- function(d) d[, k := match(paste(study, population, day, region), unique(paste(grp$study, grp$population, grp$day, grp$region)))]
wt(o(grp)[order(k, g)][, k := NULL][], a[["out_groups"]])
wt(o(gen)[order(k, g, gene)][, k := NULL][], a[["out_genes"]])
