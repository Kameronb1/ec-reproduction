# Figure S6: expression of the G3 genes and of three serotonin genes (TPH1, DDC, SLC18A1) in the published cell types of
# two human fetal pancreas studies, Olaniru 2023 (specimens, 12 to 20 weeks post-conception) and Yu 2021b (one embryo per
# age, 8 to 19 weeks), relative to random gene sets of matched expression.
#
# Profiles: raw-count pseudobulk sums per specimen or embryo and cell type, kept with at least 10 cells. In Yu 2021b the
# four endocrine progenitor labels (EP1 to EP4) are also pooled per embryo into one "Endocrine progenitors" profile.
# Gene sets: the genes of each set present in both studies (G3: 74 of 84). Score of a profile: mean log2(1 + CPM) over the
# set's genes minus the mean of the same over 500 random gene sets, each gene replaced by a gene from the same of 25 bins
# of dataset-wide average expression (average = mean of the per-cell-type means of log2(1 + CPM)); random genes exclude
# every gene-group gene. Per cell type: the mean score over specimens or embryos, and the 2.5 and 97.5% points of the
# random sets' mean scores. The pooled progenitor profiles are scored with the random sets of the per-label profiles.
#
# Inputs:  groups (gene groups), ol_counts, ol_profiles, yu_counts, yu_profiles (pseudobulk sums and profile tables)
# Outputs: points (one row per drawn specimen or embryo), rows (one row per drawn cell type)
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(data.table) })
N_NULL <- 500; MIN_CELLS <- 10; N_BINS <- 25

GG <- as.data.table(read.delim(gzfile(A[["groups"]]), check.names = FALSE, stringsAsFactors = FALSE))
gl <- function(k) GG[group == k]$gene
SER <- c("TPH1", "DDC", "SLC18A1")
FOUR <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
groups_full <- list(G1 = gl("G1_EC_GENES_SC_SIC_LACK"), G2 = gl("G2_EC_GENES_SC_SIC_HAVE"), G3 = gl("G3_SC_SIC_OWN_PROGRAM"),
                    serotonin = SER, G2_core = setdiff(gl("G2_EC_GENES_SC_SIC_HAVE"), union(SER, FOUR)))
stopifnot(identical(lengths(groups_full)[1:3], c(G1 = 28L, G2 = 34L, G3 = 84L)), length(groups_full$G2_core) == 29)
# older symbols of four current symbols in the Yu 2021b gene table; never drawn into random sets
ALIAS <- c(PAK5 = "PAK7", GALNT17 = "WBSCR17", RETREG1 = "FAM134B", STUM = "C1orf95")
EXCLUDE <- union(unlist(groups_full), ALIAS)

read_pb <- function(counts, profiles) {
  d <- read.delim(gzfile(counts, encoding = "UTF-8"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = character(),
                  colClasses = c(gene = "character"))
  meta <- as.data.table(read.delim(profiles, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8", na.strings = character(),
                                   colClasses = c(specimen = "character", label = "character")))
  M <- as.matrix(d[, -1]); storage.mode(M) <- "double"; rownames(M) <- d$gene
  stopifnot(identical(colnames(M), meta$profile), !anyDuplicated(d$gene))
  list(M = M, meta = meta)
}
keep <- function(pb, ok) { m <- pb$meta[ok]; m[, library_size := colSums(pb$M[, ok, drop = FALSE])]
  list(M = pb$M[, ok, drop = FALSE], meta = m[, .(profile, cells, specimen, label, age, library_size)]) }
OL <- read_pb(A[["ol_counts"]], A[["ol_profiles"]]); YU <- read_pb(A[["yu_counts"]], A[["yu_profiles"]])
groups <- lapply(groups_full, function(x) intersect(x, intersect(rownames(YU$M), rownames(OL$M))))
stopifnot(lengths(groups)[["G3"]] == 74, lengths(groups)[["serotonin"]] == 3)

make_bins <- function(ref) {
  expressed <- names(ref)[ref > 0]
  b <- cut(rank(ref[expressed], ties.method = "first"), N_BINS, labels = FALSE); names(b) <- expressed; b
}
draw_null <- function(bins, nrep, groups) {
  pool_genes <- setdiff(names(bins), EXCLUDE)
  pool <- split(pool_genes, bins[pool_genes])
  lapply(names(groups), function(g) {
    gg <- intersect(groups[[g]], names(bins))
    m <- replicate(nrep, vapply(gg, function(x) { p <- pool[[as.character(bins[x])]]; p[sample.int(length(p), 1)] }, ""))
    list(genes = gg, null = matrix(m, nrow = length(gg)))
  }) |> setNames(names(groups))
}
log_cpm <- function(M) log2(t(t(M) / colSums(M)) * 1e6 + 1)
null_for <- function(pb, groups) {                       # random sets from the profiles' dataset-wide average expression
  cpm <- log_cpm(pb$M); meta <- pb$meta
  ref <- rowMeans(sapply(unique(meta$label), function(p) rowMeans(cpm[, meta$label == p, drop = FALSE])))
  draw_null(make_bins(ref), N_NULL, groups)
}
score_with <- function(pb, nul) {
  cpm <- log_cpm(pb$M)
  res <- rbindlist(lapply(names(nul), function(g) {
    obs <- colMeans(cpm[nul[[g]]$genes, , drop = FALSE])
    nm <- apply(nul[[g]]$null, 2, function(s) colMeans(cpm[s, , drop = FALSE]))
    if (is.null(dim(nm))) nm <- matrix(nm, nrow = 1)
    mu <- rowMeans(nm)
    data.table(profile = colnames(cpm), group = g, delta = obs - mu, nd = split(nm - mu, row(nm)))
  }))
  merge(copy(pb$meta), res, by = "profile")
}
summ_label <- function(tab) tab[, {
  nn <- colMeans(do.call(rbind, nd))
  .(mean_delta = mean(delta), null_lo = quantile(nn, .025), null_hi = quantile(nn, .975), specimens = .N)
}, by = .(group, label)]

# Olaniru 2023
pb_ol <- keep(OL, OL$meta$cells >= MIN_CELLS)
set.seed(20261115)
tab_ol <- score_with(pb_ol, null_for(pb_ol, groups)); sl_ol <- summ_label(tab_ol)

# Yu 2021b, and its endocrine progenitors pooled per embryo
pb_yu <- keep(YU, YU$meta$cells >= MIN_CELLS)
set.seed(20261115)
nul_yu <- null_for(pb_yu, groups)
tab_yu <- score_with(pb_yu, nul_yu); sl_yu <- summ_label(tab_yu)
ep <- YU$meta[label %in% c("EP1", "EP2", "EP3", "EP4")]
EPS <- sort(unique(ep$specimen))
M_ep <- sapply(EPS, function(s) rowSums(YU$M[, ep[specimen == s]$profile, drop = FALSE]))
colnames(M_ep) <- paste0(EPS, "|Endocrine progenitors")
meta_ep <- ep[, .(cells = sum(cells), age = age[1]), by = specimen][match(EPS, specimen)]
meta_ep[, `:=`(profile = colnames(M_ep), label = "Endocrine progenitors")]
pb_ep <- keep(list(M = M_ep, meta = meta_ep), meta_ep$cells >= MIN_CELLS)
tab_ep <- score_with(pb_ep, nul_yu); sl_ep <- summ_label(tab_ep)

# drawn rows
COLS <- c(G3 = "Separate SC-SIC program (G3)", serotonin = "Serotonin genes")
OL_ROWS <- c(`endocrine progenitors` = "Endocrine progenitors", beta = "Beta", alpha = "Alpha", delta = "Delta", Schwann = "Schwann",
             mesenchymal = "Mesenchymal", acinar = "Acinar", ductal = "Ductal", endothelial = "Endothelial", immune = "Immune")
YU_ROWS <- setNames(c("Trunk", "Tip", "Endocrine progenitors", "Beta", "Alpha", "Delta", "Neuron", "Mesenchymal", "Acinar", "Ductal",
                      "Endothelial", "Immune"),
                    c("Trunk", "Tip", "Endocrine progenitors", "β", "α", "δ", "Neuron", "Mesenchyme", "Acinar cell", "Duct",
                      "Endothelium", "Immune cell"))
ENDO_ROWS <- c("Trunk", "Tip", "Endocrine progenitors", "Beta", "Alpha", "Delta")
pts <- rbind(
  tab_ol[group %in% names(COLS) & label %in% names(OL_ROWS), .(study = "ol", group, label, specimen, cells, delta)],
  tab_yu[group %in% names(COLS) & label %in% names(YU_ROWS), .(study = "yu", group, label, specimen, cells, delta)],
  tab_ep[group %in% names(COLS), .(study = "yu", group, label, specimen, cells, delta)])
bands <- rbind(
  sl_ol[group %in% names(COLS) & label %in% names(OL_ROWS), .(study = "ol", group, label, mean_delta, null_lo, null_hi, n = specimens)],
  sl_yu[group %in% names(COLS) & label %in% names(YU_ROWS), .(study = "yu", group, label, mean_delta, null_lo, null_hi, n = specimens)],
  sl_ep[group %in% names(COLS), .(study = "yu", group, label, mean_delta, null_lo, null_hi, n = specimens)])
tag <- function(d) {
  d <- copy(d)
  d[, row := ifelse(study == "ol", OL_ROWS[label], YU_ROWS[label])]
  d[, kind := fifelse(row %in% ENDO_ROWS, "endo", "other")]
  d[, col := factor(COLS[group], levels = COLS)]
  d
}
pts <- tag(pts); bands <- tag(bands)
stopifnot(nrow(pts) == 378L, nrow(bands) == 44L)
fwrite(bands[order(study, group, row)], A[["rows"]], sep = "\t")
fwrite(pts[order(study, group, row, specimen)], A[["points"]], sep = "\t")
