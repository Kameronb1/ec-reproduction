# Figure 5A: RNA levels of 34 candidate regulators (transcription factors) in each cell type, in four groups by
# association (EC / SC-SIC, endocrine, islet, intestinal).
#
# From raw-count pseudobulk sums of cells with published labels (one profile per sample and cell type), every profile above
# 569,504 reads is binomially thinned to that depth (20 draws, seed 20260925, one seeded stream per cohort), and each
# regulator's log2(1 + CPM) is averaged over the draws (tf_expression_by_sample.tsv). Intestinal EC = EC (published
# labels) of four intestinal studies (the Elmentaite 2020 pediatric Crohn's disease donors excluded); SC-SIC, SC-beta and
# SC-alpha = the six SC-islet studies; adult islet beta and alpha = two studies, adult donors only (beta_major and
# beta_minor pooled per donor; Bandesh 2026 Proliferating_Alpha pooled with Alpha per donor, as Figure 1C-D counts it;
# donor HPAP034, aged 13, left out). All islet profiles are thinned from one seeded stream, so each islet profile's draws
# depend on the other islet profiles in the stream, the dividing alpha cells included. An intestinal study with no sample
# reaching 20 EC cells (Elmentaite 2020) enters as one pooled profile of all its EC cells, thinned the same way
# (tf_expression_pooled_sparse_gut.tsv). Panel 5A table: samples with >= 20 cells, averaged within a study, studies
# averaged within a column; plus, per regulator, in how many of the six SC-islet studies SC-SIC carries more RNA than
# SC-beta (f5_5a_regulator_rna.tsv). McCauley 2023 HIO graft EC (GSE214852; EC = the published Endocrine_1 cluster) is
# computed the same way, one profile per organoid, mean of the two organoids (f5_5a_regulator_rna__hio_mccauley_rows.tsv).
# Descriptive: RNA levels say nothing about regulator activity.
#
# Inputs (name=path): <cohort>_sums / <cohort>_groups for cohort sc, gut_primary, islet, hio
#                     (pseudobulk sums: gene x profile; groups: one row per profile)
# Outputs: by_sample, pooled, panel, hio_rows
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
rt <- function(p) read.delim(if (grepl("gz$", p)) gzfile(p) else p, check.names = FALSE, stringsAsFactors = FALSE)
wt <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA")

# Regulators whose RNA is reported in the per-sample tables: five primary sets, two further sets (below) and the panel's
# TFs (SETS). Every gene is thinned, so the length of this list changes no value.
SETS_PLAN <- list(`EC selectors` = c("LMX1A", "FEV", "NKX2-2", "RFX3", "PAX4", "ASCL1", "GLIS3"),
                  `Intestinal and regional` = c("CDX1", "CDX2", "SATB2", "HOXA9", "HOXA10", "HOXB13", "PITX2", "GATA4", "GATA5", "GATA6"),
                  Pancreatic = c("PDX1", "NKX6-1", "MAFB", "ISL1", "ARX", "LMX1B"),
                  Progenitor = c("NEUROG3", "HES6", "SOX4", "INSM1", "NEUROD1"),
                  `SC-SIC associated` = c("NEUROD4", "MYT1L", "DACH2"))
SETS_EXTRA <- list(`Signal response` = c("SMAD1", "SMAD5", "SMAD9", "VDR", "NR1H4", "AHR"),
                   `Replicated programs` = c("BACH2", "ETS2", "XBP1", "SOX11", "NFE2L1"))
tfs <- unique(unlist(c(SETS_PLAN, SETS_EXTRA)))
# display groups of panel 5A: four groups by association (EC / SC-SIC, endocrine, islet, intestinal), in
# the order shown. HOXA9, HOXA10 and GATA5 stay out of the panel (they stay in the thinned set,
# which every table but the panel keeps)
SETS <- list(`EC / SC-SIC association` = c("LMX1A", "FEV", "LMX1B", "NEUROD4", "MNX1"),
             `Endocrine association` = c("NEUROG3", "INSM1", "NEUROD1", "NKX2-2", "PAX4", "RFX6", "ZNF800", "PAX6", "ARX", "ISL1",
                                         "HES6", "SOX4", "ASCL1"),
             `Islet association` = c("PDX1", "NKX6-1", "MAFA", "MAFB", "RFX3", "GLIS3", "MYT1L", "DACH2", "PLAGL1"),
             `Intestinal association` = c("GATA6", "CDX1", "CDX2", "GATA4", "SATB2", "HOXB13", "PITX2"))
TFS <- unlist(SETS, use.names = FALSE); stopifnot(length(TFS) == 34, !anyDuplicated(TFS))
tfs <- unique(c(tfs, TFS))
SC_STUDIES <- c("GSE151117", "GSE167880", "GSE185036", "GSE199636", "GSE202497", "GSE226346")
ADULT_ISLET <- c("GSE148073", "GSE221156"); ADULT_EXCLUDED <- "HPAP034"
BETA_LABELS <- c("beta", "beta_major", "beta_minor", "Beta"); ALPHA_LABELS <- c("alpha", "Alpha", "Proliferating_Alpha")
CROHNS_FAMILY <- "SF_ELMENTAITE2020_DEVGUT"; CROHNS_UNITS <- c("T017", "T019", "T176", "T189", "T197", "T202", "T203")
NAME <- c(SF_HICKEY2023_ADULT_EEC = "Hickey 2023", SF_ELMENTAITE2021_GUT_ATLAS = "Elmentaite 2021", SF_ELMENTAITE2020_DEVGUT = "Elmentaite 2020",
          SF_EGOZI2021_DEVELOPMENTAL_EEC = "Egozi 2021", SF_GSE214852_HIO_GRAFT = "McCauley 2023")
POP_LABEL <- c(SC_SIC = "SC-SIC", SC_BETA = "SC-beta", SC_ALPHA = "SC-alpha")
DEPTH <- 569504; REPS <- 20L; SEED <- 20260925L

# ---- per-sample regulator RNA (thinned) --------------------------------------------------------------------------
thin_levels <- function(m, info, cohort) {
  lib <- colSums(m); idx <- match(tfs, rownames(m)); valid <- which(!is.na(idx))
  sums <- sumsq <- matrix(0, length(tfs), ncol(m), dimnames = list(tfs, colnames(m)))
  set.seed(SEED); for (rep in seq_len(REPS)) {
    # all genes are drawn so the thinned library total (the CPM denominator) is exact
    tm <- sapply(seq_len(ncol(m)), function(j) if (lib[j] >= DEPTH) rbinom(nrow(m), m[, j], DEPTH / lib[j]) else m[, j]); den <- colSums(tm)
    lv <- matrix(NA_real_, length(tfs), ncol(m)); lv[valid, ] <- log2(1 + t(t(tm[idx[valid], , drop = FALSE]) / den) * 1e6)
    sums <- sums + lv; sumsq <- sumsq + lv^2
  }
  meanlv <- sums / REPS; mcse <- sqrt(pmax((sumsq - REPS * meanlv^2) / (REPS - 1), 0) / REPS)
  do.call(rbind, lapply(seq_len(nrow(info)), function(j) data.frame(cohort = cohort, family = info$family_id[j], unit = info$unit[j], population = info$class[j],
    cells = info$n_cells[j], library = info$lib[j], tf = tfs, mean_log2cpm = meanlv[, j], mcse = mcse[, j], min20 = info$n_cells[j] >= 20, below_depth = lib[j] < DEPTH)))
}
L <- list()
for (prefix in c("sc", "gut_primary", "islet")) {
  info <- rt(a[[paste0(prefix, "_groups")]]); x <- rt(a[[paste0(prefix, "_sums")]])
  if (prefix == "gut_primary") {
    ex <- info$family_id == CROHNS_FAMILY & info$unit %in% CROHNS_UNITS
    # the excluded profiles are all below the thinning depth, so no other profile's draws move
    stopifnot(any(ex), all(info$lib[ex] < DEPTH)); info <- info[!ex, ]
  }
  stopifnot(!anyDuplicated(x$gene), !anyDuplicated(info$group), all(info$group %in% names(x)))
  m <- as.matrix(x[, info$group]); rownames(m) <- x$gene; stopifnot(all(is.finite(m)), all(m >= 0), all(m == floor(m)), all(colSums(m) == info$lib))
  if (prefix == "islet") {   # adults only; published labels mapped to beta/alpha; beta_major + beta_minor pooled per donor
    info <- info[info$family_id %in% ADULT_ISLET & !(info$unit %in% ADULT_EXCLUDED), ]
    info$class <- ifelse(info$class %in% BETA_LABELS, "beta", ifelse(info$class %in% ALPHA_LABELS, "alpha", NA)); info <- info[!is.na(info$class), ]
    key <- paste(info$family_id, info$unit, info$class, sep = "|"); m <- sapply(split(info$group, key), function(g) rowSums(m[, g, drop = FALSE]))
    info <- do.call(rbind, lapply(split(info, key), function(z) data.frame(group = paste(z$family_id[1], z$unit[1], z$class[1], sep = "|"), family_id = z$family_id[1],
      unit = z$unit[1], class = z$class[1], n_cells = sum(z$n_cells), lib = sum(z$lib))))
    m <- m[, info$group, drop = FALSE]; stopifnot(all(colSums(m) == info$lib))
    na <- tapply(info$n_cells[info$class == "alpha"], info$family_id[info$class == "alpha"], sum)   # adult alpha cells per study
    stopifnot(na[["GSE221156"]] == 25935, na[["GSE148073"]] == 3463)
  }
  if (prefix == "sc") stopifnot(setequal(unique(info$family_id), SC_STUDIES))
  keep <- if (prefix == "sc") info$class %in% c("SC_SIC", "SC_BETA", "SC_ALPHA") else if (prefix == "gut_primary") info$class == "POSITIVE" else
    info$class %in% c("alpha", "beta")
  info <- info[keep, ]; m <- m[, info$group, drop = FALSE]
  L[[prefix]] <- thin_levels(m, info, prefix)
}
lv <- do.call(rbind, L); wt(lv, a[["by_sample"]])

# ---- gut studies with no >= 20-cell EC sample: one pooled EC profile --------------------------------------------
info <- rt(a[["gut_primary_groups"]]); x <- rt(a[["gut_primary_sums"]])
info <- info[info$class == "POSITIVE", ]; info <- info[!(info$family_id == CROHNS_FAMILY & info$unit %in% CROHNS_UNITS), ]
sparse <- names(which(tapply(info$n_cells, info$family_id, max) < 20)); stopifnot(length(sparse) > 0)
out <- list()
for (f in sparse) {
  g <- info[info$family_id == f, ]; m <- rowSums(as.matrix(x[, g$group, drop = FALSE])); names(m) <- x$gene; lib <- sum(m); stopifnot(lib == sum(g$lib), all(tfs %in% names(m)))
  set.seed(SEED); s <- ss <- setNames(numeric(length(tfs)), tfs)
  for (r in seq_len(REPS)) { t_ <- if (lib >= DEPTH) rbinom(length(m), m, DEPTH / lib) else m; names(t_) <- names(m); l2 <- log2(1 + t_[tfs] / sum(t_) * 1e6); s <- s + l2; ss <- ss + l2^2 }
  mu <- s / REPS; mcse <- sqrt(pmax((ss - REPS * mu^2) / (REPS - 1), 0) / REPS)
  out[[f]] <- data.frame(cohort = "gut_primary", family = f, unit = sprintf("pooled_%d_samples", nrow(g)), population = "POSITIVE", cells = sum(g$n_cells), library = lib,
                         tf = tfs, mean_log2cpm = mu, mcse = mcse, min20 = sum(g$n_cells) >= 20, below_depth = lib < DEPTH)
}
pz <- do.call(rbind, out); stopifnot(unique(pz$cells[pz$family == CROHNS_FAMILY]) == 45); wt(pz, a[["pooled"]])

# ---- panel 5A table ---------------------------------------------------------------------------------------------
x <- read.delim(a[["by_sample"]], stringsAsFactors = FALSE); x$min20 <- as.logical(x$min20)
x <- x[x$min20 & x$tf %in% TFS, ]
pz <- read.delim(a[["pooled"]], stringsAsFactors = FALSE); pz$min20 <- as.logical(pz$min20); pz <- pz[pz$tf %in% TFS, ]
stopifnot(!any(pz$family %in% x$family[x$cohort == "gut_primary"])); x <- rbind(x, pz[, names(x)])
GUT_ALL <- c("SF_HICKEY2023_ADULT_EEC", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_ELMENTAITE2020_DEVGUT", "SF_EGOZI2021_DEVELOPMENTAL_EEC")
GUT <- intersect(GUT_ALL, unique(x$family))
x$column <- NA_character_
x$column[x$cohort == "gut_primary"] <- paste0(NAME[x$family[x$cohort == "gut_primary"]], " gut EC")
x$column[x$cohort == "sc"] <- POP_LABEL[x$population[x$cohort == "sc"]]
x$column[x$cohort == "islet"] <- paste("adult islet", x$population[x$cohort == "islet"])
x <- x[!is.na(x$column), ]
st <- aggregate(mean_log2cpm ~ column + family + tf, x, mean)                       # samples -> study
D <- aggregate(mean_log2cpm ~ column + tf, st, mean)                                 # studies -> column
nstud <- aggregate(family ~ column, unique(st[, c("column", "family")]), length); names(nstud)[2] <- "studies"
D <- merge(D, nstud)
COLS <- c(paste0(NAME[GUT], " gut EC"), POP_LABEL, "adult islet beta", "adult islet alpha")
stopifnot(setequal(unique(D$column), COLS))
D$set <- factor(rep(names(SETS), lengths(SETS))[match(D$tf, TFS)], names(SETS))
sic <- st[st$column == "SC-SIC", ]; bet <- st[st$column == "SC-beta", ]
mm <- merge(sic, bet, by = c("family", "tf"), suffixes = c("_sic", "_beta")); mm$up <- mm$mean_log2cpm_sic > mm$mean_log2cpm_beta
K <- aggregate(up ~ tf, mm, sum); K$n <- aggregate(up ~ tf, mm, length)$up; stopifnot(all(K$n == length(SC_STUDIES)))
K$set <- factor(rep(names(SETS), lengths(SETS))[match(K$tf, TFS)], names(SETS))
wt(rbind(data.frame(kind = "level", column = D$column, set = as.character(D$set), tf = D$tf, value = D$mean_log2cpm, studies = D$studies),
         data.frame(kind = "sc_sic_above_sc_beta_studies", column = "SC-SIC vs SC-beta", set = as.character(K$set), tf = K$tf, value = K$up, studies = K$n)), a[["panel"]])

# ---- McCauley 2023 HIO graft EC ---------------------------------------------------------------------------------
wlf <- function(x, p) { con <- file(p, open = "wb"); write.table(x, con, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA", eol = "\n"); close(con) }
SET_OF <- setNames(rep(names(SETS), lengths(SETS)), TFS)
info <- rt(a[["hio_groups"]]); hx <- rt(a[["hio_sums"]]); stopifnot(!anyDuplicated(hx$gene), all(info$group %in% names(hx)))
hm <- as.matrix(hx[, info$group]); rownames(hm) <- hx$gene; stopifnot(all(hm == floor(hm)), all(colSums(hm) == info$lib))
info <- info[info$class %in% c("EC", "OTHER_EEC"), ]; hm <- hm[, info$group, drop = FALSE]
stopifnot(identical(sort(info$unit[info$class == "EC"]), c("wt1", "wt2")), sum(info$n_cells[info$class == "EC"]) == 322L, sum(info$n_cells[info$class == "OTHER_EEC"]) == 276L)
lv_h <- thin_levels(hm, info, "hio_graft")
z <- lv_h[lv_h$population == "EC" & lv_h$min20 & lv_h$tf %in% TFS, ]
h5a <- aggregate(mean_log2cpm ~ tf, aggregate(mean_log2cpm ~ family + tf, z, mean), mean)
A <- data.frame(kind = "level", column = paste0(NAME[["SF_GSE214852_HIO_GRAFT"]], " HIO graft EC"), set = unname(SET_OF[h5a$tf]), tf = h5a$tf, value = h5a$mean_log2cpm, studies = 1L)
A <- A[match(TFS, A$tf), ]; stopifnot(!anyNA(A$value)); wlf(A, a[["hio_rows"]])
