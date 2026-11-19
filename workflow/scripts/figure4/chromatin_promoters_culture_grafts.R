# Figure 4D, 4E and Supplementary Figure S8C: chromatin accessibility of the three gene groups (G1 = EC genes low in
# SC-SIC, G2 = EC genes SC-SIC have, G3 = the SC-SIC own program) in Augsornworawat 2023 (GSE199636) single-nucleus ATAC
# data.
#
#   4D  culture time within SC-SIC, SC-beta and SC-alpha (long-culture labels): weeks 2, 3 and 4, each group's median
#       promoter level minus the median of every gene with a promoter peak, in the same sample; each cell type thinned to
#       its own shallowest week. A second table keeps every group gene's own value (its level minus the same all-genes
#       median). Then months 6 and 12 and the median of the three grafts on the same scale (out_culture_ext, drawn).
#   4E, S8C  grafts: graft EC (txp_1-3) minus in vitro EC at week 2, for distal elements (4E; peaks within 50 kb of a TSS
#       of the gene that lie in no annotated promoter) and for promoters (S8C), each group's median gene minus the
#       all-genes median, against 10,000 random gene sets matched on the in vitro level; and whether the promoter and
#       distal calls agree. The in vitro side is the three week-2 samples Figure 4 compares the same grafts with
#       (sc_islet_1, sc_islet_2, s6_week2; graft labels): the graft and the three samples are thinned to their smallest
#       total (in that order, one seed), and each gene's in vitro value is the equal-weight mean of the three samples'
#       log2(1 + CPM).
#
# Inputs (name=path):
#   groups   gene group table (gene_groups_definitions.tsv.gz; columns group, gene)
#   gtf      GENCODE v44 / 10x GRCh38-2024-A genes.gtf.gz (transcript TSS)
#   mtx      GSE199636 peak x pseudobulk counts (peaks_by_sample_label.mtx)
#   peaks    its peak ids (peak_ids.txt, chr:start-end)
#   pbgroups its pseudobulk columns (sample_label_groups.tsv)
# Outputs: out_culture (4D weeks 2-4), out_culture_genes (4D weeks 2-4, per gene), out_culture_ext and
#   out_culture_ext_genes (4D as drawn: weeks 2-4, months 6 and 12, graft median), out_graft (SC-SIC),
#   out_graft_pop (4E, S8C: SC-SIC, SC-beta, SC-alpha), out_agreement
args <- commandArgs(trailingOnly = TRUE); a <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
suppressPackageStartupMessages({library(digest); library(Matrix); library(GenomicRanges); library(data.table)})
PROM <- 2000L; DIST <- 50000L; SEED <- 20260929L; NPERM <- 10000L
# every thinning is repeated REPS times from its seed's stream and each gene's value averaged over the draws before group
# medians, quartiles and random-set statistics; a draw-free profile gives the same value every draw
REPS <- 20L
mean_draws <- function(f) Reduce(`+`, lapply(seq_len(REPS), function(i) f())) / REPS
wt <- function(x, p) fwrite(x, p, sep = "\t", quote = TRUE, na = "NA")
fgz <- function(p) as.data.table(read.delim(gzfile(p), check.names = FALSE, stringsAsFactors = FALSE))
lcpm <- function(x, depth) log2(1 + x / depth * 1e6)
seed_of <- function(key, seed = SEED) seed + strtoi(substr(digest(key, algo = "xxhash32", serialize = FALSE), 1, 7), 16L)
thin <- function(x, d) if (sum(x) > d) rbinom(length(x), x, d / sum(x)) else x
GG <- fgz(a[["groups"]])[, .(group, gene)]
GS3 <- c(G1_EC_GENES_SC_SIC_LACK = "G1", G2_EC_GENES_SC_SIC_HAVE = "G2", G3_SC_SIC_OWN_PROGRAM = "G3"); GG[, g := GS3[group]]
peaks_of <- function(ids) { p <- do.call(rbind, strsplit(sub(":", "-", ids), "-")); GRanges(p[, 1], IRanges(as.integer(p[, 2]), as.integer(p[, 3]))) }

# ---- Augsornworawat 2023 matrix, transcript TSS and gene windows --------------------------------------------------------
gl <- readLines(gzfile(a[["gtf"]])); f <- strsplit(gl[grepl("\ttranscript\t", gl, fixed = TRUE)], "\t", fixed = TRUE); rm(gl)
TX <- data.table(chrom = vapply(f, `[`, "", 1), start = as.integer(vapply(f, `[`, "", 4)), end = as.integer(vapply(f, `[`, "", 5)), strand = vapply(f, `[`, "", 7),
                 gene = sub('.*gene_name "([^"]+)".*', "\\1", vapply(f, `[`, "", 9))); rm(f)
TX[, tss := ifelse(strand == "+", start, end)]; TX <- unique(TX[, .(chrom, tss, gene)])
ANYPROM <- resize(GRanges(TX$chrom, IRanges(TX$tss, width = 1)), 2 * PROM + 1L, fix = "center")
G <- fread(a[["pbgroups"]]); M <- as(readMM(a[["mtx"]]), "CsparseMatrix"); ids <- readLines(a[["peaks"]])
stopifnot(nrow(M) == length(ids), ncol(M) == nrow(G), G$group_id == paste0("g", seq_len(nrow(G))))
SGR <- peaks_of(ids); anyprom_s <- overlapsAny(SGR, ANYPROM)
tg <- GRanges(TX$chrom, IRanges(TX$tss, width = 1))
po <- findOverlaps(SGR, resize(tg, 2 * PROM + 1L, fix = "center")); do <- findOverlaps(SGR, resize(tg, 2 * DIST + 1L, fix = "center"))
WIN <- list(promoter = unique(data.table(gene = TX$gene[subjectHits(po)], peak = queryHits(po))),
            distal = unique(data.table(gene = TX$gene[subjectHits(do)], peak = queryHits(do)))[!anyprom_s[peak]])
colu <- function(layer, s, lab = "EC") which(G$scope == "primary" & G$author_layer == layer & G$sample_key == s & G$analysis_label == lab)

# ---- 4D: culture time within SC-SIC, SC-beta and SC-alpha, weeks 2-4 ----------------------------------------------------
CUL3 <- c("s6_week2", "s6_week3", "s6_week4"); POPS <- c(EC = "SC-SIC", Beta = "SC-beta", Alpha = "SC-alpha")
colp <- function(s, lab) which(G$scope == "primary" & G$author_layer == "author_label__longculture_wk234_6m_1yr" & G$sample_key == s & G$analysis_label == lab)
CB <- unlist(lapply(names(POPS), function(lab) {
  X <- lapply(setNames(CUL3, CUL3), function(s) rowSums(M[, colp(s, lab), drop = FALSE]))
  dP <- min(vapply(X, sum, 0))
  lapply(CUL3, function(s) {
    set.seed(seed_of(paste("within_pop", lab, s)))
    lv <- WIN$promoter[, .(n = sum(X[[s]][peak])), by = gene][, n := NULL]
    lv[, level := mean_draws(function() { x <- thin(X[[s]], dP); WIN$promoter[, .(n = sum(x[peak])), by = gene][lv, on = "gene", lcpm(n, dP)] })]; allm <- median(lv$level)
    list(groups = rbindlist(lapply(c("G1", "G2", "G3"), function(gg) { v <- lv[gene %in% GG[g == gg, gene], level]
           data.table(population = POPS[[lab]], sample = s, g = gg, genes = length(v), nuclei = sum(G$cells[colp(s, lab)]), events_full = sum(X[[s]]),
                      common_depth = dP, all_genes = nrow(lv), group_median = median(v), all_genes_median = allm, rel = median(v) - allm) })),
         # each group gene's own level minus the all-genes median (the per-gene values behind the boxes of each week)
         genes = merge(lv[, .(gene, level)], GG[, .(gene, g)], by = "gene")[
           , .(population = POPS[[lab]], sample = s, gene, g, level, all_genes_median = allm, rel = level - allm)])
  })
}), recursive = FALSE)
wt(rbindlist(lapply(CB, `[[`, "groups")), a[["out_culture"]])
fwrite(rbindlist(lapply(CB, `[[`, "genes")), a[["out_culture_genes"]], sep = "\t")   # unquoted, as in the published per-gene table

# ---- 4D extended: months 6 and 12 of culture and the median of the three grafts, on the weeks 2-4 scale ----------------
# Weeks 2-4 are the rows above. Months 6 and 12 (long-culture labels, as weeks 2-4) and grafts txp_1-3 (graft
# labels, in which EB nuclei count as SC-SIC) are thinned to the same cell type's weeks 2-4 common depth, or kept whole
# when shallower (depth = the depth used). Drawn: weeks, months and the graft median (median of the three grafts' group
# values; per gene, the median of its three graft values); the single grafts are kept for reference.
EXT <- list(c("s6_month6", "author_label__longculture_wk234_6m_1yr", "long-culture"), c("s6_month12", "author_label__longculture_wk234_6m_1yr", "long-culture"),
            c("txp_1", "author_label__graft_txp", "graft"), c("txp_2", "author_label__graft_txp", "graft"), c("txp_3", "author_label__graft_txp", "graft"))
cole <- function(layer, s, lab) which(G$scope == "primary" & G$author_layer == layer & G$sample_key == s & G$analysis_label == lab)
CE <- unlist(lapply(names(POPS), function(lab) {
  dP <- min(vapply(CUL3, function(s) sum(M[, colp(s, lab), drop = FALSE]), 0))
  lapply(EXT, function(e) {
    s <- e[1]; i <- cole(e[2], s, lab); stopifnot(length(i) == 1L); xf <- M[, i]; dd <- min(dP, sum(xf))
    set.seed(seed_of(paste("within_pop", lab, s)))
    lv <- WIN$promoter[, .(n = sum(xf[peak])), by = gene][, n := NULL]
    lv[, level := mean_draws(function() { x <- thin(xf, dP); WIN$promoter[, .(n = sum(x[peak])), by = gene][lv, on = "gene", lcpm(n, dd)] })]; allm <- median(lv$level)
    list(groups = rbindlist(lapply(c("G1", "G2", "G3"), function(gg) { v <- lv[gene %in% GG[g == gg, gene], level]
           data.table(population = POPS[[lab]], sample = s, g = gg, genes = length(v), nuclei = G$cells[i], events_full = sum(xf),
                      common_depth = dd, all_genes = nrow(lv), group_median = median(v), all_genes_median = allm, rel = median(v) - allm, labels = e[3]) })),
         genes = merge(lv[, .(gene, level)], GG[, .(gene, g)], by = "gene")[, .(population = POPS[[lab]], sample = s, gene, g, level, all_genes_median = allm, rel = level - allm)])
  })
}), recursive = FALSE)
TXP <- c("txp_1", "txp_2", "txp_3")
eg <- rbind(rbindlist(lapply(CB, `[[`, "groups"))[, labels := "long-culture"], rbindlist(lapply(CE, `[[`, "groups")))
en <- rbind(rbindlist(lapply(CB, `[[`, "genes")), rbindlist(lapply(CE, `[[`, "genes")))
gm <- eg[sample %in% TXP, .(sample = "graft_median", rel = median(rel), graft_mean_rel = mean(rel), labels = "graft"), by = .(population, g)]
eg <- rbind(eg, gm, fill = TRUE); eg[, drawn := !sample %in% TXP]
en <- rbind(en, en[sample %in% TXP, .(sample = "graft_median", rel = median(rel)), by = .(population, gene, g)], fill = TRUE); en[, drawn := !sample %in% TXP]
wt(eg, a[["out_culture_ext"]]); fwrite(en, a[["out_culture_ext_genes"]], sep = "\t", na = "NA")

# ---- 4E, S8C: grafts, promoter and distal, relative to the typical gene, with the matched random-gene-set null --------------
BASE <- c("sc_islet_1", "sc_islet_2", "s6_week2"); GRAFTS <- c("txp_1", "txp_2", "txp_3")
# the same comparison within each cell type (graft labels "EC", "Beta", "Alpha"); SC-SIC's seed keys have no prefix,
# SC-beta's and SC-alpha's start with their label
graft_rel <- function(lab) {
  pf <- if (lab == "EC") "" else paste0(lab, " ")
  for (s in c(BASE, GRAFTS)) stopifnot(length(colu("author_label__graft_txp", s, lab)) == 1L)
  BW <- lapply(setNames(BASE, BASE), function(s) rowSums(M[, colu("author_label__graft_txp", s, lab), drop = FALSE])); K <- length(BW)
  DD <- rbindlist(lapply(c("distal", "promoter"), function(win) {
    WA <- WIN[[win]]
    lv <- lapply(BW, function(w) WA[, .(n = sum(w[peak])), by = gene][, level := lcpm(n, sum(w))])
    for (k in seq_len(K)) stopifnot(identical(lv[[k]]$gene, lv[[1]]$gene))
    base <- data.table(gene = lv[[1]]$gene, level = Reduce(`+`, lapply(lv, `[[`, "level")) / K)   # in vitro level, full depth
    rbindlist(lapply(GRAFTS, function(tx) {
      a <- rowSums(M[, colu("author_label__graft_txp", tx, lab), drop = FALSE]); b <- BW
      set.seed(seed_of(if (win == "distal") paste0(pf, paste("graft", tx)) else paste0(pf, paste("graft", tx, win)))); dpt <- min(sum(a), vapply(b, sum, 0))
      nbc <- paste0("nb", seq_len(K))
      gsum <- function(a, b) WA[, c(list(na = sum(a[peak])), setNames(lapply(b, function(x) sum(x[peak])), nbc)), by = gene]
      gs <- gsum(a, b); gs <- gs[gs$na + Reduce(`+`, gs[, ..nbc]) > 0, .(gene)]   # genes with any accessible event before thinning
      a0 <- a; b0 <- b
      gs[, diff := mean_draws(function() {
        a <- a0; b <- b0
        if (sum(a) > dpt) a <- rbinom(length(a), a, dpt / sum(a))
        for (k in seq_len(K)) if (sum(b[[k]]) > dpt) b[[k]] <- rbinom(length(b[[k]]), b[[k]], dpt / sum(b[[k]]))
        t <- gsum(a, b)[gs, on = "gene"]
        t[, lcpm(na, dpt) - Reduce(`+`, lapply(.SD, lcpm, depth = dpt)) / K, .SDcols = nbc] })]
      gs <- merge(gs, base[, .(gene, base_level = level)], by = "gene"); gs[, bin := cut(base_level, unique(quantile(base_level, 0:10 / 10)), include.lowest = TRUE)]
      allm <- median(gs$diff)
      rbindlist(lapply(c("G1", "G2", "G3"), function(gg) {
        set_genes <- intersect(GG[g == gg, gene], gs$gene); obs <- median(gs[gene %in% set_genes, diff]); tab <- table(gs[gene %in% set_genes, bin])
        set.seed(seed_of(if (win == "distal") paste0(pf, paste("null", tx, gg)) else paste0(pf, paste("null", tx, gg, win)))); bins <- names(tab)[tab > 0]
        pools <- lapply(bins, function(bn) gs[bin == bn & !gene %in% set_genes, diff]); ns <- as.integer(tab[bins])
        null <- replicate(NPERM, median(unlist(Map(function(c0, k) c0[sample.int(length(c0), k)], pools, ns))))
        data.table(window = win, graft = tx, g = gg, genes = length(set_genes), all_genes = nrow(gs), depth = dpt, all_genes_median_diff = allm,
                   group_rel = obs - allm, null_q025_rel = quantile(null, 0.025) - allm, null_q975_rel = quantile(null, 0.975) - allm, null_median_rel = median(null) - allm,
                   perm_p_two_sided = (1 + sum(abs(null - median(null)) >= abs(obs - median(null)))) / (NPERM + 1))
      }))
    }))
  }))
  DD[, outside := fifelse(group_rel > null_q975_rel, "above", fifelse(group_rel < null_q025_rel, "below", "inside"))]
  DD[, `:=`(graft_nuclei = vapply(graft, function(s) sum(G$cells[colu("author_label__graft_txp", s, lab)]), 0),
            week2_nuclei = paste(vapply(BASE, function(s) sum(G$cells[colu("author_label__graft_txp", s, lab)]), 0), collapse = ";"))]
  DD
}
POPN <- c(EC = "SC-SIC", Beta = "SC-beta", Alpha = "SC-alpha")   # graft-layer analysis labels and their names
DP <- rbindlist(lapply(names(POPN), function(lab) cbind(population = POPN[[lab]], graft_rel(lab))))
DD <- DP[population == "SC-SIC"][, c("population", "graft_nuclei", "week2_nuclei") := NULL]
wt(DD, a[["out_graft"]]); wt(DP, a[["out_graft_pop"]])
wt(dcast(DD, graft + g ~ window, value.var = "outside")[, same := distal == promoter], a[["out_agreement"]])
