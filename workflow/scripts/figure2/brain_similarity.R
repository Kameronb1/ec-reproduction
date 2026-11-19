# Figure 2B: correlation of each sample with primary intestinal EC, adult islet beta cells and adult brain serotonin
# neurons.
#
# What it computes. Per-sample pseudobulk counts of primary intestinal EC (published labels, the seven Elmentaite 2020
# Crohn's disease donors removed), HIO graft EC, SC-SIC, SC-beta, SC-alpha and adult islet beta (beta_major + beta_minor
# summed per donor; adults only, the donor listed as excluded in islet_roster, HPAP034, removed), together with per-donor
# pseudobulk counts of adult brain serotonin neurons (Siletti 2023 cluster 397) and dopamine neurons (cluster 395).
# Units with fewer than 10 cells are dropped, and so is every unit of a group named in drop=. Only the units that remain
# enter the steps below. Genes shared by all sources, protein-coding on chr1-22/X with a genomic span; log2 CPM (edgeR
# TMM over all units, prior count 1); each unit's gene-span trend (loess on log10 span) subtracted, because brain nuclei
# carry many intronic reads; expressed genes kept (mean >= 1 or max >= 4). Pearson r of each unit with the serotonin-
# neuron mean, the primary intestinal EC mean and the adult islet beta mean (and the dopamine-neuron mean when those
# units are kept), on the 2,000 most variable genes; a brain, intestinal EC or adult beta unit is compared with the mean
# of the other units of its group. The single-nucleus SC-islet study (GSE199636) is kept as separate "(nuclei study)"
# groups. Figure 2B uses drop="Brain dopamine neurons": only the populations it draws enter the normalization.
#
# Inputs (name=path): gut_sums gut_groups sc_sums sc_groups hio_sums hio_groups islet_sums islet_groups (per-sample
#   pseudobulk sums and group tables), gut_roster (JSON with the excluded Crohn's disease donors), islet_roster (JSON,
#   adult islet datasets and excluded donors), brain_sums brain_groups (per-donor brain pseudobulk), drop (comma-separated
#   group names left out, may be empty).
# Outputs (name=path): similarity (r_serotonin, r_gut_ec, r_beta per unit, 3 decimals), similarity_unrounded (the same,
#   unrounded).
#
#   Rscript brain_similarity.R gut_sums=... ... drop="Brain dopamine neurons" similarity=... similarity_unrounded=...
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(org.Hs.eg.db); library(jsonlite) })
FLOOR <- 10L; NVAR <- 2000L
DROP <- if ("drop" %in% names(A) && nzchar(A[["drop"]])) strsplit(A[["drop"]], ",", fixed = TRUE)[[1]] else character()

# ---- 1. brain pseudobulk: per donor per cluster (adult serotonin and dopamine neurons) ----
rd <- function(f) { m <- as.data.table(read.delim(gzfile(f), check.names = FALSE)); g <- m$gene; m <- as.matrix(m[, -1]); rownames(m) <- g; m }
bg <- fread(A[["brain_groups"]])[study == "Siletti 2023" & cluster %in% c(397L, 395L)]
units_b <- sort(bg$group); SB <- rd(A[["brain_sums"]])[, units_b]
ncell_b <- bg$n_cells[match(units_b, bg$group)]
stopifnot(length(units_b) == 6L, sum(ncell_b) == 2314L)   # 3 donors x 2 clusters, 2,314 nuclei

# ---- 2. comparator pseudobulk ----
gs <- rd(A[["gut_sums"]]); gg <- fread(A[["gut_groups"]])
ss <- rd(A[["sc_sums"]]); sg <- fread(A[["sc_groups"]])
hs <- rd(A[["hio_sums"]]); hg <- fread(A[["hio_groups"]])
xs <- rd(A[["islet_sums"]]); xg <- fread(A[["islet_groups"]])
# adult islet beta per donor: beta_major + beta_minor summed (Fasolino 2022), adults only (islet_roster's excluded donors removed)
IA <- fromJSON(A[["islet_roster"]])$adult_islet
xg <- xg[class %in% c("Beta", "beta_major", "beta_minor") & family_id %in% IA$datasets & !(unit %in% names(IA$excluded_units))]
xg[, donor := paste(family_id, unit, "Beta", sep = "|")]
xs <- sapply(split(xg$group, xg$donor), function(g) rowSums(xs[, g, drop = FALSE]))
xg <- xg[, .(n_cells = sum(n_cells)), by = .(group = donor, family_id, unit)][, class := "Beta"]
stopifnot(identical(sort(colnames(xs)), sort(xg$group)))
crohns <- fromJSON(A[["gut_roster"]])$excluded_units
gg <- gg[class == "POSITIVE" & !(family_id == "SF_ELMENTAITE2020_DEVGUT" & unit %in% crohns)]
NUC <- "GSE199636"
meta <- rbind(
  gg[, .(group, grp = "Primary gut EC", study = family_id, n_cells, src = "gut")],
  hg[class == "EC", .(group, grp = "HIO EC", study = family_id, n_cells, src = "hio")],
  sg[class == "SC_SIC", .(group, grp = fifelse(family_id == NUC, "SC-SIC (nuclei study)", "SC-SIC"), study = family_id, n_cells, src = "sc")],
  sg[class == "SC_BETA", .(group, grp = fifelse(family_id == NUC, "SC-beta (nuclei study)", "SC-beta"), study = family_id, n_cells, src = "sc")],
  sg[class == "SC_ALPHA", .(group, grp = fifelse(family_id == NUC, "SC-alpha (nuclei study)", "SC-alpha"), study = family_id, n_cells, src = "sc")],
  xg[class == "Beta", .(group, grp = "Adult islet beta", study = family_id, n_cells, src = "islet")],
  data.table(group = units_b, grp = ifelse(grepl("\\|SER$", units_b), "Brain serotonin neurons", "Brain dopamine neurons"), study = "Siletti 2023", n_cells = ncell_b, src = "brain"))
stopifnot(all(DROP %in% meta$grp), !any(DROP %in% c("Brain serotonin neurons", "Primary gut EC", "Adult islet beta")))
meta <- meta[n_cells >= FLOOR & !grp %in% DROP]
# the seven Elmentaite 2020 Crohn's disease donors must be absent
stopifnot(!any(sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", meta[study == "SF_ELMENTAITE2020_DEVGUT", group]) %in% crohns))
srcm <- list(gut = gs, hio = hs, sc = ss, islet = xs, brain = SB)
genes <- Reduce(intersect, lapply(srcm, rownames))
counts <- do.call(cbind, lapply(seq_len(nrow(meta)), function(k) srcm[[meta$src[k]]][genes, meta$group[k]]))
colnames(counts) <- meta$group; stopifnot(all(counts >= 0), all(counts == round(counts)))
print(meta[, .(units = .N, cells = sum(n_cells)), by = grp])

# ---- 3. genes, log2 CPM, gene-span trend removed ----
y <- normLibSizes(DGEList(counts)); lc <- cpm(y, log = TRUE, prior.count = 1)
sym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = rownames(lc), keytype = "SYMBOL", columns = c("ENTREZID", "GENETYPE")))
sym <- as.data.table(sym)[!is.na(ENTREZID)][, .SD[1], by = SYMBOL]
loc <- merge(as.data.table(toTable(org.Hs.egCHRLOC)), as.data.table(toTable(org.Hs.egCHRLOCEND)), by = c("gene_id", "Chromosome"), allow.cartesian = TRUE)
loc <- loc[Chromosome %in% c(1:22, "X", "Y"), .(span = max(abs(end_location)) - min(abs(start_location)), chrom = Chromosome[1]), by = gene_id]
ann <- merge(sym, loc, by.x = "ENTREZID", by.y = "gene_id"); ann <- ann[GENETYPE == "protein-coding" & chrom %in% c(1:22, "X") & span > 0]
kept <- intersect(rownames(lc), ann$SYMBOL); span <- log10(ann$span[match(kept, ann$SYMBOL)]); lc <- lc[kept, ]
gm <- rowMeans(lc); fitg <- gm >= 1
trend <- sapply(seq_len(ncol(lc)), function(j) { r <- lc[, j] - gm
  f <- loess(r ~ s, data = data.frame(r = r[fitg], s = span[fitg]), span = .5, degree = 1)
  predict(f, data.frame(s = pmin(pmax(span, min(span[fitg])), max(span[fitg])))) })
lc <- lc - trend; expressed <- rowMeans(lc) >= 1 | apply(lc, 1, max) >= 4
lc <- lc[expressed, ]
cat("units", nrow(meta), "; genes: shared", length(genes), "-> protein-coding with span", length(kept), "-> expressed", nrow(lc), "\n")

# ---- 4. similarity to the neuron means ----
v <- apply(lc, 1, var); hv <- names(sort(v, decreasing = TRUE))[seq_len(min(NVAR, length(v)))]
cen <- function(g) rowMeans(lc[hv, meta$grp == g, drop = FALSE])
HAS_DA <- any(meta$grp == "Brain dopamine neurons")
ser_c <- cen("Brain serotonin neurons")
U <- data.table(group = meta$group, grp = meta$grp, study = meta$study, r_serotonin = apply(lc[hv, ], 2, cor, ser_c))
if (HAS_DA) { da_c <- cen("Brain dopamine neurons"); U[, r_dopamine := apply(lc[hv, ], 2, cor, da_c)] }
# a brain unit is not compared with a mean that contains itself: leave-one-donor-out for brain units
for (k in which(meta$src == "brain")) { same <- meta$grp == meta$grp[k] & seq_len(nrow(meta)) != k
  m <- rowMeans(lc[hv, same, drop = FALSE]); col <- if (meta$grp[k] == "Brain serotonin neurons") "r_serotonin" else "r_dopamine"
  set(U, k, col, cor(lc[hv, k], m)) }
simA <- copy(U); simA[, r_serotonin := round(r_serotonin, 3)]; if (HAS_DA) simA[, r_dopamine := round(r_dopamine, 3)]

# ---- 4 (continued). similarity to the primary intestinal EC mean (an intestinal EC unit against the mean of the
# other intestinal EC units) and to the adult islet beta mean (an adult beta donor against the mean of the other donors)
gut <- meta$grp == "Primary gut EC"
U[, r_gut_ec := apply(lc[hv, ], 2, cor, rowMeans(lc[hv, gut, drop = FALSE]))]
for (k in which(gut)) set(U, k, "r_gut_ec", cor(lc[hv, k], rowMeans(lc[hv, gut & seq_len(nrow(meta)) != k, drop = FALSE])))
beta <- meta$grp == "Adult islet beta"
U[, r_beta := apply(lc[hv, ], 2, cor, rowMeans(lc[hv, beta, drop = FALSE]))]
for (k in which(beta)) set(U, k, "r_beta", cor(lc[hv, k], rowMeans(lc[hv, beta & seq_len(nrow(meta)) != k, drop = FALSE])))
simA[, `:=`(r_gut_ec = round(U$r_gut_ec, 3), r_beta = round(U$r_beta, 3))]
fwrite(simA, A[["similarity"]], sep = "\t")
fwrite(U, A[["similarity_unrounded"]], sep = "\t")
cat("brain similarity tables written\n")
