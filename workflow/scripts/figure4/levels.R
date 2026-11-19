# Figure 4B (right, the four markers' expression change) and 4C (gap to primary intestinal EC on the Figure 2D gene
# groups), from pseudobulk levels of all cells (no downsampling, "unthinned").
#
# Profiles: one per sample x population (SC-SIC, SC-beta, SC-alpha) with >= 20 cells, raw UMI sums from published
# labels; Augsornworawat 2023 from its own recount (nuclei). Level = log2(1 + CPM) over the gene universe G (genes tested
# in the donor-level Figure 2C model that are present in every tree and in the intestinal reference); the B bars use G
# plus CFAP92. Graft rows (A-HUES8 ... P-TXP) use SoupX-corrected sums for both ends; the three in vitro culture rows use
# uncorrected sums. B-0-M3's later end (Balboa 2022 3-month grafts) is one pooled profile per population: its three
# graft samples' sums added (>= 20 cells pooled), SRR13805905 uncorrected (SoupX could not estimate its contamination).
# Row ends and samples come from the row table of ec_calls.R (a pooled end is its samples joined with "+").
#   B bars  (later mean - earlier mean), averaged over the four classifier markers (TPH1, SLC18A1, FEV, LMX1A)
#   C       per profile, median over a gene group of (level - primary intestinal EC reference); reference = mean of three
#           intestinal studies, (mean(Elmentaite 2021, Elmentaite 2020) + Hickey 2023 + Egozi 2021) / 3, pooled EC cells
#           with published labels; row end = mean over its samples. Drawn: the Figure 2D groups (28/34/84 genes) with
#           the seven Elmentaite 2020 Crohn's donor units left out of the reference (from/to). For comparison the table
#           also carries the same values on uncorrected counts (*_unthinned_uncorrected).
#   C bars  (E_end_spread) the spread of the group's genes at the later end: each gene's (level - reference) averaged
#           over the later samples, then the 25th, 50th and 75th percentiles over the group's genes (R quantile type 7).
#           The pale bar spans q25 to q75; the arrow end (mean over samples of per-sample medians) need not equal q50.
#
# Inputs (name=path): tree1..tree6 _profiles/_genes/_counts, gut=<intestinal reference analysis.rds>, cfap92=<CFAP92 alias row>,
#   universe=<genes tested in the Figure 2C model>, nuc_profiles, nuc_pseudobulk, soupx_profiles, soupx_pseudobulk,
#   rows=<row_samples.tsv>, groups=<Figure 2D gene groups>, roster=<Crohn's units JSON>
# Outputs: bars=C_bars.tsv  gaps=E_group_gaps.tsv  end_spread=E_end_spread.tsv
args <- commandArgs(trailingOnly = TRUE)
A <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages({ library(Matrix); library(jsonlite) })
rt <- function(p) read.delim(if (grepl("gz$", p)) gzfile(p) else p, stringsAsFactors = FALSE, check.names = FALSE)
wtsv <- function(x, f) { con <- file(f, "wb"); write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n"); close(con) }
# The per-gene E gaps are written as text (write.table, 15 significant digits) and read back before they are averaged,
# so the spread values are identical to the published ones.
via_text <- function(x) { f <- tempfile(fileext = ".tsv"); on.exit(unlink(f)); wtsv(x, f); rt(f) }
lev <- function(M) log2(1 + t(t(M) / colSums(M)) * 1e6)

MIN_CELLS <- 20L; MATURE <- c("SC_SIC", "SC_BETA", "SC_ALPHA", "SC_DELTA", "SC_EPSILON", "SC_GAMMA"); DRAWN <- c("SC_SIC", "SC_BETA", "SC_ALPHA")
FAMILIES <- c(GCA2021 = "SF_ELMENTAITE2021_GUT_ATLAS", HICKEY2023 = "SF_HICKEY2023_ADULT_EEC",
              EGOZI2021 = "SF_EGOZI2021_DEVELOPMENTAL_EEC", ELM2020 = "SF_ELMENTAITE2020_DEVGUT")
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
# pooled ends: the members' profiles added, per population
pool_profiles <- function(M, INFO, pooled) {
  members <- strsplit(pooled, "+", fixed = TRUE)[[1]]; I0 <- INFO; s <- sub("\\|.*", "", I0$id)   # members from the unpooled table
  for (p in unique(I0$population[s %in% members])) { ids <- I0$id[s %in% members & I0$population == p]
    M <- cbind(M, matrix(rowSums(M[, ids, drop = FALSE]), dimnames = list(rownames(M), paste(pooled, p, sep = "|"))))
    INFO <- rbind(INFO, data.frame(id = paste(pooled, p, sep = "|"), population = p, cells = sum(I0$cells[match(ids, I0$id)]))) }
  list(M = M, INFO = INFO) }

# ---- primary intestinal EC pools, without the seven Elmentaite 2020 Crohn's donor units ----
a <- readRDS(A[["gut"]]); md <- a$metadata; cnt <- a$counts; rownames(cnt) <- a$genes
xr <- readRDS(A[["cfap92"]]); stopifnot(identical(colnames(xr$counts), colnames(cnt)))
ec_all <- md$author_group %in% c("EXPLICIT_EC", "EXPLICIT_AUTHOR_EC") & md$family_id %in% FAMILIES
ROSTER <- fromJSON(A[["roster"]]); stopifnot(identical(ROSTER$family_id, FAMILIES[["ELM2020"]]), length(ROSTER$excluded_units) == 7L)
crohn <- md$family_id == FAMILIES[["ELM2020"]] & md$biological_unit_id %in% ROSTER$excluded_units
stopifnot(sum(ec_all & crohn) == 10L)
pools <- function(ec) { PRI <- sapply(names(FAMILIES), function(s) Matrix::rowSums(cnt[, ec & md$family_id == FAMILIES[[s]], drop = FALSE]))
  rbind(PRI, CFAP92 = sapply(names(FAMILIES), function(s) sum(xr$counts[1, ec & md$family_id == FAMILIES[[s]]]))) }
PRI_NC <- pools(ec_all & !crohn)    # drawn reference
rm(a, cnt, xr); invisible(gc())

U <- readLines(A[["universe"]])
drawn <- INFO$id[INFO$population %in% DRAWN & INFO$cells >= MIN_CELLS]
G <- intersect(U, intersect(g10, rownames(PRI_NC))); d10 <- intersect(drawn, colnames(M10))
INFOP <- rt(A[["nuc_profiles"]]); pb <- rt(A[["nuc_pseudobulk"]])
MP <- rowsum(as.matrix(pb[, INFOP$id]), pb$gene); rm(pb); stopifnot(all(colSums(MP) == INFOP$library_matched))
drawnP <- INFOP$id[INFOP$population %in% DRAWN & INFOP$cells >= MIN_CELLS]
GP <- intersect(U, intersect(rownames(MP), rownames(PRI_NC)))
G2 <- c(G, "CFAP92"); GP2 <- c(GP, "CFAP92")

# ---- SoupX-corrected profiles of the graft rows' captures ----
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

# ---- rows ----
RC <- rt(A[["rows"]])[, c("id", "ref_samples", "test_samples")]
GRAFT <- c("A-HUES8", "A-WS4", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-TXP"); CULT <- c("B-w0-w6", "P-w2-m6", "P-w2-m12")
stopifnot(setequal(RC$id, c(GRAFT, CULT)))
for (pooled in grep("+", unique(c(RC$ref_samples, RC$test_samples)), fixed = TRUE, value = TRUE)) {
  x <- pool_profiles(M10, INFO, pooled); M10 <- x$M; INFO <- x$INFO; pool_soupx(pooled) }
drawn <- INFO$id[INFO$population %in% DRAWN & INFO$cells >= MIN_CELLS]; d10 <- intersect(drawn, colnames(M10))
spl <- function(x) strsplit(x, ";", fixed = TRUE)[[1]]
REF <- setNames(lapply(RC$ref_samples, spl), RC$id); TEST <- setNames(lapply(RC$test_samples, spl), RC$id)
is_p <- function(id) grepl("^P", id)
graft_ids <- unique(unlist(lapply(GRAFT, function(id) c(REF[[id]], TEST[[id]]))))

# ---- level tables: uncorrected (U) and graft profiles SoupX-corrected (S) ----
build <- function(M, rows, prof, P) {
  U <- cbind(M[rows, prof], P[rows, ]); S <- U
  hit <- intersect(prof, spr$id[sub("\\|.*", "", spr$id) %in% graft_ids])
  S[, hit] <- SXC[rows, hit]
  list(U = lev(U), S = lev(S))
}
tables <- function(P) list(wc_G = build(M10, G, d10, P), p_G = build(MP, GP, drawnP, P), wc_G2 = build(M10, G2, d10, P), p_G2 = build(MP, GP2, drawnP, P))
L <- tables(PRI_NC)
soupx <- function(id) id %in% GRAFT
tab_for <- function(L, id, univ) { t <- L[[paste0(if (is_p(id)) "p_" else "wc_", univ)]]; if (soupx(id)) t$S else t$U }
tab_unc <- function(L, id, univ) L[[paste0(if (is_p(id)) "p_" else "wc_", univ)]]$U
cl <- function(s, p) paste(s, p, sep = "|")
chg <- function(T, id, p, genes) { r <- intersect(cl(REF[[id]], p), colnames(T)); t <- intersect(cl(TEST[[id]], p), colnames(T))
  if (!length(r) || !length(t)) return(rep(NA_real_, length(genes))); rowMeans(T[genes, t, drop = FALSE]) - rowMeans(T[genes, r, drop = FALSE]) }

# ---- B bars: the four markers' expression change ----
MK <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
CB <- do.call(rbind, lapply(c(GRAFT, CULT), function(id) do.call(rbind, lapply(c("SC_SIC", "SC_BETA"), function(p) {
  data.frame(id = id, population = p, counts = if (soupx(id)) "soupx" else "uncorrected",
    four_marker = mean(chg(tab_for(L, id, "G2"), id, p, MK)), four_marker_unthinned_uncorrected = mean(chg(tab_unc(L, id, "G2"), id, p, MK))) }))))
wtsv(CB, A[["bars"]])

# ---- C: gaps to the primary intestinal EC reference ----
GRPS <- c(G1 = "G1_EC_GENES_SC_SIC_LACK", G2 = "G2_EC_GENES_SC_SIC_HAVE", G3 = "G3_SC_SIC_OWN_PROGRAM")
GG <- rt(A[["groups"]]); stopifnot(identical(as.vector(table(GG$group)[GRPS]), c(28L, 34L, 84L)))
prim3 <- function(T) (rowMeans(T[, c("GCA2021", "ELM2020")]) + T[, "HICKEY2023"] + T[, "EGOZI2021"]) / 3
gap_arm <- function(T, id, p, which, G, GG) { s <- if (which == "ref") REF[[id]] else TEST[[id]]; cols <- intersect(cl(s, p), colnames(T))
  if (!length(cols)) return(NA_real_); g <- intersect(GG$gene[GG$group == GRPS[[G]]], rownames(T)); PR <- prim3(T)
  mean(vapply(cols, function(c) median(T[g, c] - PR[g]), 0)) }
E <- do.call(rbind, lapply(names(GRPS), function(G) do.call(rbind, lapply(c(GRAFT, CULT), function(id) do.call(rbind, lapply(DRAWN, function(p)
  data.frame(group = G, id = id, population = p, counts = if (soupx(id)) "soupx" else "uncorrected",
    from = gap_arm(tab_for(L, id, "G"), id, p, "ref", G, GG), to = gap_arm(tab_for(L, id, "G"), id, p, "test", G, GG),
    from_unthinned_uncorrected = gap_arm(tab_unc(L, id, "G"), id, p, "ref", G, GG), to_unthinned_uncorrected = gap_arm(tab_unc(L, id, "G"), id, p, "test", G, GG))))))))
wtsv(E, A[["gaps"]])

# ---- C bars: spread of the group's genes at the later end (drawn groups and reference, as from/to) ----
gap_genes <- function(T, id, p, G) { cols <- intersect(cl(TEST[[id]], p), colnames(T)); if (!length(cols)) return(NULL)
  g <- intersect(GG$gene[GG$group == GRPS[[G]]], rownames(T)); PR <- prim3(T)
  do.call(rbind, lapply(cols, function(c) data.frame(group = G, id = id, population = p, sample = sub("\\|.*$", "", c), gene = g, gap = unname(T[g, c] - PR[g])))) }
EG <- via_text(do.call(rbind, lapply(names(GRPS), function(G) do.call(rbind, lapply(c(GRAFT, CULT), function(id) do.call(rbind, lapply(DRAWN, function(p)
  gap_genes(tab_for(L, id, "G"), id, p, G))))))))
EE <- aggregate(gap ~ group + id + population + gene, EG, mean)                       # each gene's mean over the later samples
key <- paste(EE$group, EE$id, EE$population)
ES <- do.call(rbind, lapply(split(EE, factor(key, levels = sort(unique(key), method = "radix"))), function(x)
  data.frame(group = x$group[1], id = x$id[1], population = x$population[1], genes = nrow(x),
             q25 = unname(quantile(x$gap, .25)), q50 = median(x$gap), q75 = unname(quantile(x$gap, .75)))))
stopifnot(nrow(ES) == sum(!is.na(E$to)))
wtsv(ES, A[["end_spread"]])
