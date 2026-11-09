# Figure S12A, S12B, S12C: likeness of SC-islet cells to adult brain serotonin neurons.
#
# What it computes. Per-sample pseudobulk counts of primary intestinal EC (primary study labels, Elmentaite 2020 Crohn's
# donors removed), HIO graft EC, induced EC, SC-SIC, SC-beta, SC-alpha and adult islet beta (beta_major + beta_minor summed
# per donor; adults only, the donor listed as excluded in islet_roster, HPAP034, removed), together with per-donor pseudobulk
# counts of adult brain serotonin neurons (Siletti 2023 cluster 397) and dopamine neurons (cluster 395). Units with fewer
# than 10 cells are dropped. Genes shared by all sources, protein-coding on chr1-22/X with a genomic span; log2 CPM
# (edgeR TMM over all units, prior count 1); each unit's gene-span trend (loess on log10 span) subtracted, because brain
# nuclei carry many intronic reads; expressed genes kept (mean >= 1 or max >= 4).
#   A  Pearson r of each unit with the serotonin-neuron mean, the dopamine-neuron mean, the primary intestinal EC mean and
#      the adult islet beta mean, on the 2,000 most variable genes (a brain, intestinal EC or adult beta unit is compared
#      with the mean of the other units of its group).
#   B  serotonin-versus-dopamine lean: SER genes are >= 2 log2 higher in serotonin than dopamine neurons in every donor
#      pair, DA genes the reverse; score = mean(SER genes) - mean(DA genes), genes centered on the all-unit mean.
#   C  each SER gene's mean level in SC-SIC, SC-beta, primary intestinal EC and serotonin neurons.
# The single-nucleus SC-islet study (GSE199636) is kept as separate "(nuclei study)" groups in A and B; C pools all
# ten SC-islet samples.
#
# Inputs (name=path): gut_sums gut_groups sc_sums sc_groups hio_sums hio_groups ind_sums ind_groups islet_sums
#   islet_groups (per-sample pseudobulk sums and group tables), gut_roster (JSON with the excluded Crohn's donors),
#   islet_roster (JSON, adult islet datasets and excluded donors),
#   brain_sums brain_groups (per-donor brain pseudobulk).
# Outputs (name=path): similarity (A: r_serotonin, r_dopamine, r_gut_ec, r_beta per unit), units (B: per-unit lean, with
#   A's serotonin and dopamine correlations),
#   summary (B: group medians), ser_da_genes (the SER and DA gene lists), ser_genes (C).
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(org.Hs.eg.db); library(jsonlite) })
FLOOR <- 10L; NVAR <- 2000L; T_LOG2 <- 2

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
is_ <- rd(A[["ind_sums"]]); ig <- fread(A[["ind_groups"]])
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
  ig[class == "EC", .(group, grp = "Induced EC", study = family_id, n_cells, src = "ind")],
  sg[class == "SC_SIC", .(group, grp = fifelse(family_id == NUC, "SC-SIC (nuclei study)", "SC-SIC"), study = family_id, n_cells, src = "sc")],
  sg[class == "SC_BETA", .(group, grp = fifelse(family_id == NUC, "SC-beta (nuclei study)", "SC-beta"), study = family_id, n_cells, src = "sc")],
  sg[class == "SC_ALPHA", .(group, grp = fifelse(family_id == NUC, "SC-alpha (nuclei study)", "SC-alpha"), study = family_id, n_cells, src = "sc")],
  xg[class == "Beta", .(group, grp = "Adult islet beta", study = family_id, n_cells, src = "islet")],
  data.table(group = units_b, grp = ifelse(grepl("\\|SER$", units_b), "Brain serotonin neurons", "Brain dopamine neurons"), study = "Siletti 2023", n_cells = ncell_b, src = "brain"))
meta <- meta[n_cells >= FLOOR]
srcm <- list(gut = gs, hio = hs, ind = is_, sc = ss, islet = xs, brain = SB)
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
lc <- lc - trend; expressed <- rowMeans(lc) >= 1 | apply(lc, 1, max) >= 4; lc <- lc[expressed, ]
cat("genes: shared", length(genes), "-> protein-coding with span", length(kept), "-> expressed", nrow(lc), "\n")

# ---- 4A. similarity to the two neuron means (the intestinal EC and adult beta means follow at the end) ----
v <- apply(lc, 1, var); hv <- names(sort(v, decreasing = TRUE))[seq_len(min(NVAR, length(v)))]
cen <- function(g) rowMeans(lc[hv, meta$grp == g, drop = FALSE])
ser_c <- cen("Brain serotonin neurons"); da_c <- cen("Brain dopamine neurons")
simA <- data.table(group = meta$group, grp = meta$grp, study = meta$study,
                   r_serotonin = round(apply(lc[hv, ], 2, cor, ser_c), 3), r_dopamine = round(apply(lc[hv, ], 2, cor, da_c), 3))
# a brain unit is not compared with a mean that contains itself: leave-one-donor-out for brain units
for (k in which(meta$src == "brain")) { same <- meta$grp == meta$grp[k] & seq_len(nrow(meta)) != k
  m <- rowMeans(lc[hv, same, drop = FALSE]); col <- if (meta$grp[k] == "Brain serotonin neurons") "r_serotonin" else "r_dopamine"
  set(simA, k, col, round(cor(lc[hv, k], m), 3)) }

# ---- 4B. serotonin-versus-dopamine lean ----
sb <- meta$grp == "Brain serotonin neurons"; db <- meta$grp == "Brain dopamine neurons"
d_min <- apply(lc, 1, function(z) min(outer(z[sb], z[db], "-"))); d_max <- apply(lc, 1, function(z) max(outer(z[sb], z[db], "-")))
ser_genes <- rownames(lc)[d_min >= T_LOG2]; da_genes <- rownames(lc)[d_max <= -T_LOG2]
cat("SER genes:", length(ser_genes), " DA genes:", length(da_genes), "\n")
cz <- lc - rowMeans(lc)
simB <- data.table(group = meta$group, grp = meta$grp, study = meta$study,
                   ser_mean = round(colMeans(cz[ser_genes, ]), 3), da_mean = round(colMeans(cz[da_genes, ]), 3))
simB[, ser_minus_da := ser_mean - da_mean]
res <- merge(simA, simB, by = c("group", "grp", "study"))
fwrite(res, A[["units"]], sep = "\t")
fwrite(data.table(gene = c(ser_genes, da_genes), set = rep(c("SER", "DA"), c(length(ser_genes), length(da_genes)))), A[["ser_da_genes"]], sep = "\t")
summ <- res[, .(units = .N, r_serotonin = round(median(r_serotonin), 3), r_dopamine = round(median(r_dopamine), 3),
                ser_minus_da = round(median(ser_minus_da), 3)), by = grp]
fwrite(summ, A[["summary"]], sep = "\t"); print(summ)

# ---- 4C. SER genes in SC-SIC against SC-beta and intestinal EC ----
# every sample of the group, the single-nucleus SC-islet study included (as in Figure S12D and Figure 3D)
gmean <- function(g) rowMeans(lc[, meta$grp %in% c(g, paste(g, "(nuclei study)")), drop = FALSE])
shared <- data.table(gene = ser_genes, ser_neuron = round(gmean("Brain serotonin neurons")[ser_genes], 2), sc_sic = round(gmean("SC-SIC")[ser_genes], 2),
                     sc_beta = round(gmean("SC-beta")[ser_genes], 2), gut_ec = round(gmean("Primary gut EC")[ser_genes], 2))
shared[, `:=`(sic_minus_beta = sc_sic - sc_beta, sic_minus_gut = sc_sic - gut_ec)]
fwrite(shared[order(-sic_minus_beta)], A[["ser_genes"]], sep = "\t")

# ---- 4A (continued). similarity to the primary intestinal EC mean (an intestinal EC unit against the mean of the
# other intestinal EC units) ----
gut <- meta$grp == "Primary gut EC"
simA[, r_gut_ec := round(apply(lc[hv, ], 2, cor, rowMeans(lc[hv, gut, drop = FALSE])), 3)]
for (k in which(gut)) set(simA, k, "r_gut_ec", round(cor(lc[hv, k], rowMeans(lc[hv, gut & seq_len(nrow(meta)) != k, drop = FALSE])), 3))
# ... and to the adult islet beta mean (an adult beta donor against the mean of the other donors)
beta <- meta$grp == "Adult islet beta"
simA[, r_beta := round(apply(lc[hv, ], 2, cor, rowMeans(lc[hv, beta, drop = FALSE])), 3)]
for (k in which(beta)) set(simA, k, "r_beta", round(cor(lc[hv, k], rowMeans(lc[hv, beta & seq_len(nrow(meta)) != k, drop = FALSE])), 3))
fwrite(simA, A[["similarity"]], sep = "\t")
cat("brain similarity tables written\n")
