# Figure 4F: the fetal-to-adult age score of Figure 3H (age_scores.py), for SC-SIC and SC-beta in every start and end
# sample of Figure 4C's rows, and for Figure 3H's seven intestinal rows (the ticks on 4F's axis).
#
# A profile is the pooled raw counts of its own cells, as log2(CPM + 1), with CPM over the 21,520 genes the intestinal
# pseudobulk measures. Each list gene (Figure 3H's list with the 10% detection floor: 7 adult-high, 9 fetal-high) is
# centered on the mean of the four intestinal arm profiles (Elmentaite 2020 + 2021 fetal, the mean of the two parts;
# Egozi 2021 fetal; Elmentaite 2021 adult; Hickey 2023). Score = mean centered value of the adult-high genes minus that of
# the fetal-high genes; 0 lies midway between fetal and adult intestinal EC. Every profile needs >= 20 cells. The Balboa
# 2022 3-month grafts are scored pooled (sample named by its members joined with "+"), as in Figure 4C.
#
# Inputs: tree1..tree5 _profiles/_genes/_counts (Balboa 2022, Augsornworawat 2020 HUES8, WS4corr, Du 2022 stage 6,
#   Du 2022 graft), nuc_profiles, nuc_pseudobulk (Augsornworawat 2023 recount), gut_sums, gut_groups (intestinal EEC by
#   primary study label), test_sums, test_groups (intestinal test sources by four-marker classifier result), adult_floor, fetal_floor,
#   rows (rows_called.tsv)
# Outputs: out=<age_score_by_sample.tsv>  ticks=<fig3H_ticks.tsv>  genes=<each list gene's values per sample and population>
args <- commandArgs(trailingOnly = TRUE)
A <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
rt <- function(p) read.delim(if (grepl("gz$", p)) gzfile(p) else p, stringsAsFactors = FALSE, check.names = FALSE)
wt <- function(x, f) { con <- file(f, "wb"); on.exit(close(con)); write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n") }
MIN_CELLS <- 20L
LISTS <- list(adult = scan(A[["adult_floor"]], what = "", quiet = TRUE), fetal = scan(A[["fetal_floor"]], what = "", quiet = TRUE))
GENES <- c(LISTS$adult, LISTS$fetal)

# ---- intestinal profiles: the gene universe, the per-gene center and the seven tick rows ----
pb <- function(tag) { s <- rt(A[[paste0(tag, "_sums")]]); g <- rt(A[[paste0(tag, "_groups")]])
  m <- as.matrix(s[, -1]); rownames(m) <- s$gene; stopifnot(identical(colnames(m), g$group)); list(m = m, g = g) }
GUT <- pb("gut"); TST <- pb("test")
UNIVERSE <- rownames(GUT$m); stopifnot(all(GENES %in% UNIVERSE), all(UNIVERSE %in% rownames(TST$m)))
own <- function(v, cells) { stopifnot(cells >= MIN_CELLS, all(UNIVERSE %in% names(v))); v <- v[UNIVERSE]; log2(v[GENES] / sum(v) * 1e6 + 1) }
pool <- function(p, mask) own(rowSums(p$m[, mask, drop = FALSE]), sum(p$g$n_cells[mask]))
gg <- GUT$g
arm <- function(fam, unit_ok = rep(TRUE, nrow(gg))) pool(GUT, gg$family_id == fam & gg$class == "POSITIVE" & unit_ok)
ARMS <- list(`Elmentaite 2020 + 2021 fetal` = (arm("SF_ELMENTAITE2020_DEVGUT", grepl("^[0-9]", gg$unit)) +
                                               arm("SF_ELMENTAITE2021_GUT_ATLAS", startsWith(gg$unit, "F"))) / 2,
             `Egozi 2021 fetal` = arm("SF_EGOZI2021_DEVELOPMENTAL_EEC", endsWith(gg$unit, "w")),
             `Elmentaite 2021 adult` = arm("SF_ELMENTAITE2021_GUT_ATLAS", startsWith(gg$unit, "A")),
             `Hickey 2023 adult (nuclei)` = arm("SF_HICKEY2023_ADULT_EEC"))
CENTER <- Reduce(`+`, ARMS) / length(ARMS)
score <- function(level) { x <- level - CENTER; mean(x[LISTS$adult]) - mean(x[LISTS$fetal]) }
TEST <- c(`Yu 2021 fetal *` = "SF_YU2021_ENDODERM", `Wang 2020 adult *` = "SF_GSE125970_ADULT_INTESTINE",
          `Burclaff 2022 adult *` = "SF_GSE185224_ADULT_REGIONAL_INTESTINE")
tick_levels <- c(ARMS, lapply(TEST, function(f) pool(TST, TST$g$family_id == f & TST$g$class == "CALLED")))
TICKS <- data.frame(row = names(tick_levels), kind = rep(c("intestinal arm", "intestinal test source"), c(4, 3)),
                    group = ifelse(grepl("fetal", names(tick_levels)), "fetal", "adult"), score = vapply(tick_levels, score, 0))
wt(TICKS, A[["ticks"]])

# ---- SC-SIC and SC-beta per sample: whole-cell trees and the Augsornworawat 2023 recount ----
STUDY <- c("Balboa2022", "Augsornworawat2020", "Augsornworawat2020", "Du2022", "Du2022")
PROF <- list(); INFO <- list()
for (i in seq_along(STUDY)) {
  pr <- rt(A[[sprintf("tree%d_profiles", i)]]); pr <- pr[pr$variant == "primary" & pr$population %in% c("SC_SIC", "SC_BETA"), ]
  ge <- rt(A[[sprintf("tree%d_genes", i)]]); ct <- rt(A[[sprintf("tree%d_counts", i)]])
  m <- as.matrix(ct[, pr$profile, drop = FALSE]); rownames(m) <- ge$gene[match(ct$feature_index, ge$feature_index)]
  m <- rowsum(m, rownames(m)); s <- sub(".*__", "", pr$sample)
  for (k in seq_len(nrow(pr))) PROF[[paste(s[k], pr$population[k])]] <- list(v = m[, k], cells = pr$cells[k])
  INFO[[i]] <- data.frame(study = STUDY[i], sample = s, context = pr$context)
}
NI <- rt(A[["nuc_profiles"]]); np <- rt(A[["nuc_pseudobulk"]]); MN <- rowsum(as.matrix(np[, NI$id]), np$gene)
for (s in unique(NI$sample)) for (p in c("SC_SIC", "SC_BETA")) { x <- NI$id[NI$sample == s & NI$population == p]
  if (length(x)) PROF[[paste(s, p)]] <- list(v = rowSums(MN[, x, drop = FALSE]), cells = sum(NI$cells[NI$id %in% x])) }
INFO <- unique(rbind(do.call(rbind, INFO), unique(NI[, c("study", "sample", "context")])))

RC <- rt(A[["rows"]]); RC <- RC[RC$rule == "FOUR_MARKER_UNGATED" & RC$population == "SC_SIC", ]
SAMPLES <- unique(unlist(strsplit(c(RC$ref_samples, RC$test_samples), ";")))
OUT <- list(); GN <- list()
for (s in SAMPLES) for (p in c("SC_SIC", "SC_BETA")) {
  parts <- PROF[paste(strsplit(s, "+", fixed = TRUE)[[1]], p)]; stopifnot(!any(vapply(parts, is.null, TRUE)))
  level <- own(Reduce(`+`, lapply(parts, `[[`, "v")), sum(vapply(parts, `[[`, 0, "cells")))
  x <- level - CENTER
  OUT[[length(OUT) + 1]] <- data.frame(sample = s, population = p, age_score = score(level), adult_genes = length(LISTS$adult),
                                       fetal_genes = length(LISTS$fetal), cells = sum(vapply(parts, `[[`, 0, "cells")))
  GN[[length(GN) + 1]] <- data.frame(sample = s, population = p, gene = GENES, list = rep(c("adult-high", "fetal-high"), lengths(LISTS)),
                                     log2cpm = unname(level), centered = unname(x))
}
member <- function(s) strsplit(s, "+", fixed = TRUE)[[1]][1]
AG <- do.call(rbind, OUT); k <- match(vapply(AG$sample, member, ""), INFO$sample); stopifnot(!anyNA(k))
AG <- cbind(AG[, c("sample")], INFO[k, c("study", "context")], AG[, -1]); names(AG)[1] <- "sample"
wt(AG, A[["out"]])
wt(do.call(rbind, GN), A[["genes"]])
