# Figure S1B: how many fitting cells in each fitting group end up in each four-marker class, and which
# step removes them (score cutoff, then local consistency; the 2-of-4 detection floor is kept as a
# check and removes no cell at or above the cutoff).
#
# Inputs : cells (per-cell classes written by marker_fit.R)
# Output : out (one row per fitting group, largest first, then all groups together)
#   Rscript outcome_by_group.R cells=... out=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))

CLASS_HIGH <- "HIGH_CONFIDENCE_EC"; CLASS_LOW <- "LOW_CONFIDENCE_EC"; CLASS_NONE <- "INSUFFICIENT_EC_EVIDENCE"
UNIT_OF <- c("E-MTAB-10187_intestinal" = "SF_YU2021_ENDODERM:fetal", "E-MTAB-8901_fetal" = "SF_ELMENTAITE2020_DEVGUT:fetal",
             "E-MTAB-8901_paediatric" = "SF_ELMENTAITE2020_DEVGUT:child", "GSE125970" = "SF_GSE125970_ADULT_INTESTINE",
             "GSE185224" = "SF_GSE185224_ADULT_REGIONAL_INTESTINE")
NAME_OF <- c("E-MTAB-10187_intestinal" = "Yu 2021 fetal", "E-MTAB-8901_fetal" = "Elmentaite 2020 fetal",
             "E-MTAB-8901_paediatric" = "Elmentaite 2020 child", "GSE125970" = "Wang 2020", "GSE185224" = "Burclaff 2022")

d <- read.delim(gzfile(a[["cells"]]), stringsAsFactors = FALSE)
stopifnot(nrow(d) == 1199L, setequal(unique(d$arm_id), names(UNIT_OF)))
total_high <- sum(d$presentation_class == CLASS_HIGH)
stopifnot(sum(d$score_pass & !d$floor_pass) == 0L)          # the undrawn floor removes nobody

arm_n <- table(d$arm_id)
arm_order <- names(sort(arm_n, decreasing = TRUE))
wide <- do.call(rbind, lapply(c(arm_order, "ALL_ARMS"), function(g) {
  all <- identical(g, "ALL_ARMS")
  s <- if (all) d else d[d$arm_id == g, ]
  entering <- c(nrow(s), sum(s$candidate))
  removed <- c(sum(!s$score_pass), sum(s$candidate & !s$locally_consistent))
  data.frame(scope = if (all) "all_arms" else "arm", arm_id = g,
             unit_id = if (all) NA_character_ else unname(UNIT_OF[[g]]),
             study = if (all) NA_character_ else unname(NAME_OF[[g]]),
             cells = nrow(s),
             high_confidence_ec = sum(s$presentation_class == CLASS_HIGH),
             low_confidence_ec = sum(s$presentation_class == CLASS_LOW),
             insufficient_ec_evidence = sum(s$presentation_class == CLASS_NONE),
             high_pct_of_arm = 100 * sum(s$presentation_class == CLASS_HIGH) / nrow(s),
             high_pct_of_all_high = 100 * sum(s$presentation_class == CLASS_HIGH) / total_high,
             score_cutoff_entering = entering[1], score_cutoff_removed = removed[1],
             score_cutoff_pct_removed = 100 * removed[1] / entering[1],
             local_consistency_entering = entering[2], local_consistency_removed = removed[2],
             local_consistency_pct_removed = 100 * removed[2] / entering[2],
             detection_floor_entering = as.integer(sum(s$score_pass)),
             detection_floor_removed = as.integer(sum(s$score_pass & !s$floor_pass)),
             detection_floor_drawn = FALSE,
             retained_high_confidence = sum(s$locally_consistent),
             stringsAsFactors = FALSE)
}))
write.table(wide, a[["out"]], sep = "\t", row.names = FALSE, quote = TRUE)
