# Figure S16A: the 26 measured cholesterol biosynthesis genes of Figure 4D in each graft comparison and over in vitro culture time
# alone, per gene and as the median over the genes, for SC-SIC, SC-beta and SC-alpha.
#
# Levels as Figure 4D (levels.R): per sample x population profile (>= 20 cells) of raw UMI sums from primary study labels,
# log2(1 + CPM) over a gene universe, all cells (no downsampling). Graft rows use SoupX-corrected sums for both ends
# (the pooled 3-month grafts: their three samples' sums added, SRR13805905 uncorrected, as in Figure 4D); culture rows
# use uncorrected sums. Universes: G (genes tested in the donor-level Figure 3D model that are present in every whole-cell tree
# and in the intestinal reference) for the whole-cell studies and Zhu 2023; GP (tested genes in the Augsornworawat 2023 recount
# and the intestinal reference) for the nuclei; GV (tested genes in the Veres 2019 profiles, its pseudobulk tree and the intestinal
# reference) for Veres 2019. Rows: the eleven Figure 4D comparisons (samples from the Figure 4C rows table) and three more
# in vitro steps from the comparison table: Balboa 2022 week 0 -> 3, Zhu 2023 day 21 -> 39 and Veres 2019 week 0 -> 5.
# Zhu 2023 and Veres 2019 use separately built profiles that include SC-alpha (same source files and label rule; their
# SC-SIC and SC-beta counts are checked against the pseudobulk trees).
# Per gene: ref_mean, test_mean, change (test - ref), pair_min / pair_max = smallest and largest change over every
# earlier x later sample pair; dot = >= 2 pairs, |change| >= 0.5 and every pair moving the same way; gene_column = the
# gene's rank by its mean change over the graft rows (largest first). Medians over the 26 genes per row and population.
# fill (tile color) and lab (bar label) are the drawn values, so the tables hold what the panel shows; fill uses
# the scales package (installed with ggplot2).
# Check: the eleven Figure 4D comparisons equal Figure 4D's per-gene changes and medians.
#
# Inputs (name=path): tree1..tree7 _profiles/_genes/_counts, gut=<intestinal reference analysis.rds (gene names only)>,
#   universe=<genes tested in the Figure 3D model>, nuc_profiles, nuc_pseudobulk, soupx_profiles, soupx_pseudobulk,
#   rows=<rows_called.tsv>, comparisons=<comparison definitions, whole-cell studies>, chol=<Reactome cholesterol biosynthesis genes>,
#   zhu_alpha_profiles/_counts, veres_alpha_profiles/_counts, d_genes=<D_cholesterol_per_gene.tsv>,
#   d_medians=<D_cholesterol_medians.tsv>
# Outputs: tiles=<per gene table>  medians=<median per row and population>
args <- commandArgs(trailingOnly = TRUE)
A <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages(library(Matrix))
rt <- function(p) read.delim(if (grepl("gz$", p)) gzfile(p) else p, stringsAsFactors = FALSE, check.names = FALSE)
write_lf <- function(x, path) { con <- file(path, "wb"); on.exit(close(con)); write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n") }
lev <- function(M) log2(1 + t(t(M) / colSums(M)) * 1e6)

MIN_CELLS <- 20L; MATURE <- c("SC_SIC", "SC_BETA", "SC_ALPHA", "SC_DELTA", "SC_EPSILON", "SC_GAMMA"); DRAWN <- c("SC_SIC", "SC_BETA", "SC_ALPHA")
GRAFT <- c("A-HUES8", "A-WS4", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-TXP"); CULT <- c("B-w0-w6", "P-w2-m6", "P-w2-m12")
ALL <- c(GRAFT, "B-w0-w3", "B-w0-w6", "Z-D21-D39", "V-w0-w5", "P-w2-m6", "P-w2-m12")
spl <- function(x) strsplit(x, ";", fixed = TRUE)[[1]]

# ---- whole-cell trees (universe G), the Augsornworawat 2023 recount (GP) ----
read_tree <- function(i) {
  pr <- rt(A[[sprintf("tree%d_profiles", i)]]); pr <- pr[pr$variant == "primary" & pr$population %in% MATURE, ]
  ge <- rt(A[[sprintf("tree%d_genes", i)]]); ct <- rt(A[[sprintf("tree%d_counts", i)]])
  m <- as.matrix(ct[, pr$profile, drop = FALSE]); rownames(m) <- ge$gene[match(ct$feature_index, ge$feature_index)]
  stopifnot(all(colSums(m) == pr$library_size))
  m <- rowsum(m, rownames(m)); colnames(m) <- paste(sub(".*__", "", pr$sample), pr$population, sep = "|")
  list(m = m, info = data.frame(id = colnames(m), population = pr$population, cells = pr$cells))
}
TR <- lapply(1:6, read_tree); INFO <- do.call(rbind, lapply(TR, `[[`, "info")); stopifnot(!anyDuplicated(INFO$id))
g10 <- Reduce(intersect, lapply(TR, function(x) rownames(x$m))); M10 <- do.call(cbind, lapply(TR, function(x) x$m[g10, , drop = FALSE])); rm(TR)
# A pooled end ("a+b+c" in the rows table) gets one profile per population: its samples' raw sums added, as in levels.R.
pool_profiles <- function(M, INFO, pooled) {
  members <- strsplit(pooled, "+", fixed = TRUE)[[1]]; I0 <- INFO; s <- sub("\\|.*", "", I0$id)
  for (p in unique(I0$population[s %in% members])) { ids <- I0$id[s %in% members & I0$population == p]
    M <- cbind(M, matrix(rowSums(M[, ids, drop = FALSE]), dimnames = list(rownames(M), paste(pooled, p, sep = "|"))))
    INFO <- rbind(INFO, data.frame(id = paste(pooled, p, sep = "|"), population = p, cells = sum(I0$cells[match(ids, I0$id)]))) }
  list(M = M, INFO = INFO) }
RC0 <- rt(A[["rows"]])
for (pooled in unique(grep("+", unlist(strsplit(c(RC0$ref_samples, RC0$test_samples), ";", fixed = TRUE)), fixed = TRUE, value = TRUE))) {
  x <- pool_profiles(M10, INFO, pooled); M10 <- x$M; INFO <- x$INFO }
GUT_GENES <- readRDS(A[["gut"]])$genes; invisible(gc())
U <- readLines(A[["universe"]])
drawn <- INFO$id[INFO$population %in% DRAWN & INFO$cells >= MIN_CELLS]
G <- intersect(U, intersect(g10, GUT_GENES)); d10 <- intersect(drawn, colnames(M10))
INFOP <- rt(A[["nuc_profiles"]]); pb <- rt(A[["nuc_pseudobulk"]])
MP <- rowsum(as.matrix(pb[, INFOP$id]), pb$gene); rm(pb); stopifnot(all(colSums(MP) == INFOP$library_matched))
drawnP <- INFOP$id[INFOP$population %in% DRAWN & INFOP$cells >= MIN_CELLS]
GP <- intersect(U, intersect(rownames(MP), GUT_GENES))

# ---- SoupX-corrected profiles replace the graft rows' profiles ----
SAMPLE_OF <- c(SRR10359526 = "GSM4143264", s6_week2 = "s6_week2_graft_labels")
spr <- rt(A[["soupx_profiles"]])
spr$id <- paste(ifelse(spr$capture %in% names(SAMPLE_OF), SAMPLE_OF[spr$capture], spr$capture), spr$population, sep = "|")
x <- rt(A[["soupx_pseudobulk"]]); SXC <- as.matrix(x[, spr$profile]); rownames(SXC) <- x$gene; colnames(SXC) <- spr$id; rm(x)
# a pooled graft end's SoupX profile: its captures' corrected sums added (B-0-M3: SRR13805905 enters uncorrected, as
# SoupX could not estimate that capture's contamination; its rows in the SoupX tables carry its raw counts)
pool_soupx <- function(pooled) {
  members <- strsplit(pooled, "+", fixed = TRUE)[[1]]
  for (p in unique(spr$population[spr$capture %in% members])) {
    ids <- spr$id[spr$capture %in% members & spr$population == p]
    SXC <<- cbind(SXC, matrix(rowSums(SXC[, ids, drop = FALSE]), dimnames = list(rownames(SXC), paste(pooled, p, sep = "|"))))
    spr <<- rbind(spr, transform(spr[match(ids[1], spr$id), ], profile = NA, capture = pooled, id = paste(pooled, p, sep = "|"))) }
}
for (pooled in unique(grep("+", unlist(strsplit(c(RC0$ref_samples, RC0$test_samples), ";", fixed = TRUE)), fixed = TRUE, value = TRUE))) pool_soupx(pooled)
RC <- rt(A[["rows"]]); RC <- unique(RC[RC$rule == "FOUR_MARKER_UNGATED" & RC$population == "SC_SIC", c("id", "ref_samples", "test_samples")])
stopifnot(setequal(RC$id, c(GRAFT, CULT)))
REF <- setNames(lapply(RC$ref_samples, spl), RC$id); TEST <- setNames(lapply(RC$test_samples, spl), RC$id)
graft_ids <- unique(unlist(lapply(GRAFT, function(id) c(REF[[id]], TEST[[id]]))))
build <- function(M, rows, prof) {
  U <- M[rows, prof]; S <- U
  hit <- intersect(prof, spr$id[sub("\\|.*", "", spr$id) %in% graft_ids])
  S[, hit] <- SXC[rows, hit]
  list(U = lev(U), S = lev(S))
}
L <- list(wc = build(M10, G, d10), p = build(MP, GP, drawnP))

# ---- the three extra in vitro steps ----
CM <- rt(A[["comparisons"]])
for (id in c("B-w0-w3", "Z-D21-D39", "V-w0-w5")) { k <- CM$comparison == id; stopifnot(sum(k) == 1L); REF[[id]] <- spl(CM$ref[k]); TEST[[id]] <- spl(CM$test[k]) }
stopifnot(identical(REF[["B-w0-w3"]], REF[["B-w0-w6"]]))
read_alpha <- function(k) {
  pr <- rt(A[[paste0(k, "_profiles")]]); ct <- rt(A[[paste0(k, "_counts")]]); m <- as.matrix(ct[, -1]); rownames(m) <- ct$gene
  stopifnot(!anyDuplicated(ct$gene), all(pr$cells >= MIN_CELLS))
  colnames(m) <- sub("^[^_]+__", "", colnames(m)); m }
MZ <- read_alpha("zhu_alpha"); stopifnot(all(G %in% rownames(MZ)))      # Zhu 2023 on universe G
zsb <- grep("SC_(SIC|BETA)$", colnames(MZ), value = TRUE); stopifnot(all(zsb %in% colnames(M10)), all(MZ[G, zsb] == M10[G, zsb]))
L$zhu <- list(U = lev(MZ[G, ]))
MV <- read_alpha("veres_alpha")                                          # Veres 2019 on its own universe GV
vpr <- rt(A[["tree7_profiles"]]); vge <- rt(A[["tree7_genes"]]); vct <- rt(A[["tree7_counts"]])
vpr <- vpr[vpr$variant == "primary" & vpr$population %in% c("SC_SIC", "SC_BETA"), ]
mv <- as.matrix(vct[, vpr$profile]); rownames(mv) <- vge$gene[match(vct$feature_index, vge$feature_index)]
rownames(mv)[rownames(mv) == "FAM134B"] <- "RETREG1"; mv <- rowsum(mv, rownames(mv))
colnames(mv) <- paste(sub(".*__", "", vpr$sample), vpr$population, sep = "|")
GV <- intersect(U, intersect(intersect(rownames(MV), rownames(mv)), GUT_GENES))
vsb <- grep("SC_(SIC|BETA)$", colnames(MV), value = TRUE); stopifnot(all(vsb %in% colnames(mv)), all(MV[GV, vsb] == mv[GV, vsb]))
L$veres <- list(U = lev(MV[GV, ]))

# ---- per-gene effects ----
CHT <- rt(A[["chol"]]); stopifnot(nrow(CHT) == 27L)
CH <- sort(intersect(CHT$gene, Reduce(intersect, list(G, GP, GV))), method = "radix"); stopifnot(length(CH) == 26L)
tab <- function(id) if (id == "Z-D21-D39") L$zhu$U else if (id == "V-w0-w5") L$veres$U else { t <- L[[if (grepl("^P", id)) "p" else "wc"]]; if (id %in% GRAFT) t$S else t$U }
EF <- do.call(rbind, lapply(ALL, function(id) do.call(rbind, lapply(DRAWN, function(p) {
  T <- tab(id); stopifnot(all(CH %in% rownames(T)))
  r <- intersect(paste(REF[[id]], p, sep = "|"), colnames(T)); t <- intersect(paste(TEST[[id]], p, sep = "|"), colnames(T))
  if (!length(r) || !length(t)) return(NULL)
  X <- T[CH, , drop = FALSE]; pair <- matrix(sapply(t, function(tt) sapply(r, function(rr) X[, tt] - X[, rr])), nrow = length(CH))
  data.frame(comparison = id, population = p, gene = CH,
             ref_mean = rowMeans(X[, r, drop = FALSE]), test_mean = rowMeans(X[, t, drop = FALSE]),
             change = rowMeans(X[, t, drop = FALSE]) - rowMeans(X[, r, drop = FALSE]), pair_min = apply(pair, 1, min), pair_max = apply(pair, 1, max),
             ref_n = length(r), test_n = length(t), row.names = NULL) }))))
stopifnot(nrow(EF) == length(ALL) * 3L * 26L)
MED <- aggregate(change ~ comparison + population, EF, median)

# check: the eleven Figure 4D comparisons equal Figure 4D's per-gene changes and medians
DG <- rt(A[["d_genes"]]); DM <- rt(A[["d_medians"]])
k <- match(paste(DG$comparison, DG$population, DG$gene), paste(EF$comparison, EF$population, EF$gene))
km <- match(paste(DM$comparison, DM$population), paste(MED$comparison, MED$population))
stopifnot(!anyNA(k), !anyNA(km), max(abs(EF$change[k] - DG$change)) < 1e-12, max(abs(MED$change[km] - DM$change)) < 1e-12)

# ---- S16A tables ----
EF$dot <- EF$ref_n * EF$test_n >= 2 & abs(EF$change) >= 0.5 & ((EF$pair_min > 0 & EF$change > 0) | (EF$pair_max < 0 & EF$change < 0))
EF$group <- ifelse(EF$comparison %in% GRAFT, "graft", "in_vitro")
GO <- names(sort(tapply(EF$change[EF$group == "graft"], EF$gene[EF$group == "graft"], mean), decreasing = TRUE))
EF$gene_column <- match(EF$gene, GO)
L3 <- 3; PAL <- scales::pal_div_gradient("#2A7F8E", "grey97", "#B4553F", "Lab")
EF$fill <- PAL(scales::rescale_mid(pmin(pmax(EF$change, -L3), L3), to = c(0, 1), from = c(-L3, L3), mid = 0))
ord <- order(match(EF$comparison, ALL), match(EF$population, DRAWN), EF$gene_column)
write_lf(EF[ord, c("comparison", "population", "gene", "gene_column", "change", "dot", "fill", "ref_n", "test_n", "group",
                   "ref_mean", "test_mean", "pair_min", "pair_max")], A[["tiles"]])
lab1 <- function(v) { s <- sprintf("%.1f", v); ifelse(s %in% c("-0.0", "0.0"), "0.0", ifelse(v > 0, paste0("+", s), sub("-", "\u2212", s, fixed = TRUE))) }
MED$group <- ifelse(MED$comparison %in% GRAFT, "graft", "in_vitro"); MED$lab <- lab1(MED$change)
write_lf(MED[order(match(MED$comparison, ALL), match(MED$population, DRAWN)), c("comparison", "population", "group", "change", "lab")], A[["medians"]])
