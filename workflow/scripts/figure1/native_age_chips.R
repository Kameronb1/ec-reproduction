# Figure 1A: age ranges of the native-intestine donors, one chip per life stage.
# Each native-intestine study part has a curated age record (low and high age with units, taken from the sample
# records or the paper). For every stage (fetal, neonatal, child, adult) the chip spans the smallest low and the
# largest high age across the records of that stage, converted to the stage's unit (weeks, days, years).
# Input:  ages    curated age records, one row per study part
# Output: out     chips table (stage, min, max, record giving the min and the max, chip label)
#   Rscript native_age_chips.R ages=<native_age_records.tsv> out=<native_age_chips.tsv>
args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
stopifnot(all(c("ages", "out") %in% names(arg)))

ag <- read.delim(arg[["ages"]], colClasses = "character", na.strings = character(0), encoding = "UTF-8")
ag$study <- gsub("Yu Q 2021", "Yu 2021", ag$study, fixed = TRUE)   # paper-wide study name
stopifnot(length(unique(ag$study)) == 6, all(c(ag$low_unit, ag$high_unit) %in% c("day", "week", "year")))

to_unit <- function(v, u, out) { v <- as.numeric(v); d <- v * c(day = 1, week = 7, year = 365.25)[u]
  unname(d / c(day = 1, week = 7, year = 365.25)[out]) }
STAGE <- data.frame(stage = c("fetal", "neonatal", "child", "adult"), unit = c("week", "day", "year", "year"),
                    abbr = c("wk", "d", "y", "y"), stringsAsFactors = FALSE)
stopifnot(all(ag$stage %in% STAGE$stage), all(STAGE$stage %in% ag$stage))
CHIPS <- do.call(rbind, lapply(seq_len(nrow(STAGE)), function(k) { r <- ag[ag$stage == STAGE$stage[k], ]
  lo <- min(to_unit(r$low, r$low_unit, STAGE$unit[k])); hi <- max(to_unit(r$high, r$high_unit, STAGE$unit[k]))
  data.frame(stage = STAGE$stage[k], min = lo, max = hi, min_from = paste(r$study, r$part)[which.min(to_unit(r$low, r$low_unit, STAGE$unit[k]))],
             max_from = paste(r$study, r$part)[which.max(to_unit(r$high, r$high_unit, STAGE$unit[k]))],
             label = sprintf("%s  %d–%d %s", STAGE$stage[k], round(lo), round(hi), STAGE$abbr[k]), stringsAsFactors = FALSE) }))
write.table(CHIPS, arg[["out"]], sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
