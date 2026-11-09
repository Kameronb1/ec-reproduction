# Figure 4B and Figure S14A: mature endocrine composition of every SC-islet sample, from primary study labels.
#
# For each sample, cells are counted per primary-label category (SC-SIC, SC-alpha, SC-beta, SC-delta, other hormone
# cells; proliferating alpha cells are "cycling"; everything else, including the other cells labeled as dividing, is
# "not mature endocrine"), and each mature category is given as a percentage of the sample's mature endocrine cells.
# Timing and replicate labels come from the publications (table TIMING below). Samples with fewer than 50 mature
# endocrine cells are flagged "low count". Figure 4B additionally flags the Augsornworawat 2023 month 6 and month 12
# in vitro samples ("outlier, not interpreted"); Figure S14A does not. In Figure 4B the Augsornworawat 2023 stage 6
# week 2 sample is counted on the graft-object labels, as are SC-islet 1 and 2 and Figures 4C, 4E and 4F; Figure S14A
# counts it on the long-culture labels.
#
# Inputs: per-sample profile tables of the Balboa 2022, Augsornworawat 2020 (HUES8, WS4corr) and Du 2022 trees;
# the Zhu 2023 per-cell primary study labels; the Augsornworawat 2023 recount profiles (per-sample nuclei per label) and
# its primary study label groups; counts for Veres 2019 and Siehler 2024.
# Outputs: f4b=<table for Figure 4B>  s1=<table for Figure S14A>
#
#   Rscript composition.R balboa=... hues8=... ws4=... du_s6d2=... du_10wpt=... zhu=... a23_profiles=... \
#           a23_groups=... extra=... f4b=... s1=...
args <- commandArgs(trailingOnly = TRUE)
A <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages(library(data.table))

# ---- Augsornworawat 2023 nuclei per sample and label (whole-islet culture series and grafts) ----
POPS <- c("SC_SIC", "SC_BETA", "SC_ALPHA", "SC_DELTA")
SERIES <- c("s6_week2", "s6_week3", "s6_week4", "s6_month6", "s6_month12"); GRAFTS <- c("txp_1", "txp_2", "txp_3")
INFO <- read.delim(A[["a23_profiles"]], stringsAsFactors = FALSE)
INFO <- INFO[INFO$population %in% POPS & INFO$sample %in% c(SERIES, GRAFTS), ]
AB <- do.call(rbind, lapply(c(SERIES, GRAFTS), function(s) {
  x <- INFO[INFO$sample == s, ]
  data.frame(sample = s, context = x$context[1], population = POPS, nuclei = sapply(POPS, function(p) sum(x$cells[x$population == p])))
}))

# ---- counts per sample x category ----
TREES <- c(Balboa2022 = A[["balboa"]], Augsornworawat2020 = A[["hues8"]], Augsornworawat2020 = A[["ws4"]],
           Du2022 = A[["du_s6d2"]], Du2022 = A[["du_10wpt"]])
CAT <- c("SC-SIC", "SC-alpha", "SC-beta", "SC-delta", "other hormone")
TREE_MAP <- c(SC_SIC = "SC-SIC", SC_ALPHA = "SC-alpha", SC_BETA = "SC-beta", SC_DELTA = "SC-delta", SC_EPSILON = "other hormone", SC_GAMMA = "other hormone")
rows <- list()
for (i in seq_along(TREES)) {
  p <- fread(TREES[i])[variant == "primary"]
  rows[[length(rows) + 1]] <- p[, .(study = names(TREES)[i], sample = sub(".*__", "", sample), context,
                                    category = ifelse(population %in% names(TREE_MAP), TREE_MAP[population], "not mature endocrine"), cells)]
}
z <- as.data.table(read.delim(gzfile(A[["zhu"]]), stringsAsFactors = FALSE))[author_label_status == "AUTHOR_LABELLED"]
ZM <- c(SC_EC = "SC-SIC", SC_alpha = "SC-alpha", SC_beta = "SC-beta", SC_delta = "SC-delta")
rows[[length(rows) + 1]] <- z[, .(cells = .N), by = .(sample = capture_id, context = stage_day, category = ifelse(author_anno %in% names(ZM), ZM[author_anno], "not mature endocrine"))][, study := "Zhu2023"]
e <- fread(A[["extra"]])
g <- fread(A[["a23_groups"]])[scope == "primary"]
GM <- c(EC = "SC-SIC", Alpha = "SC-alpha", Beta = "SC-beta", Delta = "SC-delta", `Proliferating Alpha` = "cycling")
AM <- c(SC_SIC = "SC-SIC", SC_ALPHA = "SC-alpha", SC_BETA = "SC-beta", SC_DELTA = "SC-delta")

# Augsornworawat 2023: SC-islet 1 and 2 carry the graft-object labels (author_label__graft_txp); the culture series and
# the grafts come from the recount profiles, with proliferating alpha nuclei from the long-culture labels as "cycling".
# week2_graft = TRUE (Figure 4B): the stage 6 week 2 sample also takes every label from the
# graft-object layer, as in Figures 4C, 4E and 4F; its long-culture proliferating-alpha row is left out (the graft object
# has no such label). Its rows keep the category order of the long-culture rows (SC-SIC, SC-beta, SC-alpha, SC-delta).
# week2_graft = FALSE (Figure S14A): week 2 on the long-culture labels.
# pool_m3 = TRUE (Figure 4B): the three Balboa 2022 3-month grafts are counted as one pooled sample, as in Figures
# 4A and 4C-F (1 to 33 SC-SIC cells each, too few to read singly); Figure S14A lists them one by one.
M3 <- c("SRR13805903", "SRR13805904", "SRR13805905")
counts <- function(week2_graft, pool_m3) {
  W2 <- "s6_week2"
  a <- as.data.table(AB)
  if (week2_graft) a <- a[sample != W2]
  r <- c(rows, list(a[, .(study = "Augsornworawat2023", sample, context, category = AM[population], cells = nuclei)]))
  graft_rows <- function(keys) g[author_layer == "author_label__graft_txp" & sample_key %in% keys,
    .(study = "Augsornworawat2023", sample = sample_key, context = ifelse(sample_key == W2, "S6_week2", sample_key),
      category = ifelse(analysis_label %in% names(GM), GM[analysis_label], "not mature endocrine"), cells)]
  r[[length(r) + 1]] <- graft_rows(c("sc_islet_1", "sc_islet_2"))
  if (week2_graft) {
    w2 <- graft_rows(W2)
    W2_ORDER <- c(unname(AM[POPS]), "not mature endocrine")             # the long-culture rows' order
    stopifnot(nrow(w2) == 5L, setequal(w2$category, W2_ORDER))
    r[[length(r) + 1]] <- w2[order(match(category, W2_ORDER))]
  }
  r[[length(r) + 1]] <- g[author_layer == "author_label__longculture_wk234_6m_1yr" & analysis_label == "Proliferating Alpha" & !(week2_graft & sample_key == W2),
                          .(study = "Augsornworawat2023", sample = sample_key, context = AB$context[match(sample_key, AB$sample)], category = "cycling", cells)]
  stopifnot(!anyNA(r[[length(r)]]$context))
  r[[length(r) + 1]] <- e[, .(study, sample, context, category, cells)]
  r <- rbindlist(r, use.names = TRUE)
  if (pool_m3) {
    stopifnot(setequal(r[study == "Balboa2022" & context == "M3_graft", unique(sample)], M3))
    r[study == "Balboa2022" & context == "M3_graft", `:=`(sample = paste(M3, collapse = "+"), context = "M3_graft_pooled")]
  }
  r[, .(cells = sum(cells)), by = .(study, sample, context, category)]
}

# ---- timing and row labels (from each publication) ----
TIMING <- rbindlist(list(
  data.table(study = "Veres2019", context = c(paste0("S5_day_", 0:7), paste0("S6_week_", 0:5)), ord = 1:14,
             when = c(paste0("S5 day ", 0:7), paste0("S6 week ", 0:5)), day = c(13:20, NA, 24, 31, 38, 45, 52)),
  data.table(study = "Augsornworawat2020", context = c("S6_in_vitro", "WS4_baseline", "M6_graft_HUES8", "M6_graft_WS4corr"), ord = 1:4,
             when = c("S6 day 14 (HUES8)", "S6 (WS4corr)", "6-mo graft (HUES8)", "6-mo graft (WS4corr)"), day = c(33, NA, NA, NA)),
  data.table(study = "Balboa2022", context = c("S5_DT", "S5_NT", "S7w0", "S7w3", "S7w6", "M1_graft", "M3_graft", "M6_graft_CTRL", "M6_graft_STZ"), ord = 1:9,
             when = c("S5 end, DT", "S5 end, NT", "S7 week 0", "S7 week 3", "S7 week 6", "1-mo graft (from S7 wk 0)", "3-mo graft (from S7 wk 0)", "6-mo graft, control (from S7 wk 3)", "6-mo graft, STZ (from S7 wk 3)"),
             day = c(16, 16, 24, 45, 66, NA, NA, NA, NA)),
  data.table(study = "Balboa2022", context = "M3_graft_pooled", ord = 7, when = "3-mo grafts, pooled (from S7 wk 0)", day = NA),
  data.table(study = "Du2022", context = c("S6D2", "10_wpt"), ord = 1:2, when = c("S6 day 2", "10-week graft"), day = c(23, NA)),
  data.table(study = "Zhu2023", context = c("D11", "D14", "D21", "D32", "D39"), ord = 1:5,
             when = c("day 11", "day 14", "day 21", "day 32", "day 39"), day = c(11, 14, 21, 32, 39)),
  data.table(study = "Augsornworawat2023", context = c("sc_islet_1", "sc_islet_2", "S6_week2", "S6_week3", "S6_week4", "S6_month6", "S6_month12", "M6_graft"), ord = 1:8,
             when = c("SC-islet 1", "SC-islet 2", "S6 week 2", "S6 week 3", "S6 week 4", "S6 month 6", "S6 month 12", "6-mo graft"),
             day = c(33, 33, 33, 40, 47, 201, 384, NA)),
  data.table(study = "Siehler2024", context = "S6_WT_D40", ord = 1, when = "S6 WT", day = 40)))
summarise <- function(AL) {
  AL <- merge(AL, TIMING, by = c("study", "context"), all.x = TRUE)
  if (anyNA(AL$ord)) stop("no timing row for: ", paste(unique(AL[is.na(ord), paste(study, context)]), collapse = ", "))
  AL[, arm := ifelse(grepl("graft|wpt", context, ignore.case = TRUE), "graft", "in dish")]
  reps <- unique(AL[, .(study, context, sample)])[order(study, context, sample)][, rep := seq_len(.N), by = .(study, context)][, nrep := .N, by = .(study, context)]
  AL <- merge(AL, reps, by = c("study", "context", "sample"))
  AL[study == "Veres2019", rep := as.integer(sub("differentiation_([0-9]).*", "\\1", sample))]
  AL[, row_lab := ifelse(nrep > 1 | study == "Veres2019", paste0(when, "  #", rep), when)]
  S <- AL[, .(mature = sum(cells[category %in% CAT]), cycling = sum(cells[category == "cycling"]), other = sum(cells[!category %in% c(CAT, "cycling")])),
          by = .(study, sample, context, row_lab, ord, rep, arm, day)]
  P <- merge(AL[category %in% CAT], S[, .(study, sample, mature)], by = c("study", "sample"))[, pct := 100 * cells / mature]
  list(S = S, P = P)
}
LOW_N <- 50
write_table <- function(P, S, path) {
  fwrite(merge(P[, .(study, sample, context, category, cells, pct)], S[, .(study, sample, row_lab, arm, day, mature, cycling, not_mature_other = other, flag)], by = c("study", "sample")),
         path, sep = "\t", eol = "\n")
}
B <- summarise(counts(week2_graft = TRUE, pool_m3 = TRUE))
S4B <- copy(B$S)[, flag := fifelse(study == "Augsornworawat2023" & context %in% c("S6_month6", "S6_month12"), "outlier, not interpreted",
                            fifelse(mature < LOW_N, "low count", fifelse(study == "Du2022" & arm == "graft", "source sample not stated", "")))]
write_table(B$P, S4B, A[["f4b"]])
S1 <- summarise(counts(week2_graft = FALSE, pool_m3 = FALSE))
SS1 <- copy(S1$S)[, flag := fifelse(mature < LOW_N, "low count", fifelse(study == "Du2022" & arm == "graft", "source sample not stated", ""))]
write_table(S1$P, SS1, A[["s1"]])
