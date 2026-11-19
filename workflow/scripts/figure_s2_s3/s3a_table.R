# Figure S3A: the values drawn. For each test study (three intestinal studies and the McCauley 2023 HIO grafts),
# the share of EC-labeled cells and of other terminal EEC classified as EC at every cutoff of the 0.02 grid and at
# the cutoff used ("curve"), and the five leave-one-fitting-group-out refits, each at its position on the full
# fit's score scale with the share of EC-labeled cells it classifies as EC ("refit").
#
# Inputs : constants (fitted_constants.json), sweep and refits (heldout.R), hio_sweep and hio_refits (hio_four_marker.R)
# Output : out
#   Rscript s3a_table.R constants=... sweep=... refits=... hio_sweep=... hio_refits=... out=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
rd <- function(p) read.delim(if (grepl("[.]gz$", p)) gzfile(p) else p, stringsAsFactors = FALSE, check.names = FALSE)

CUT <- jsonlite::fromJSON(a[["constants"]])$ec_core_cutoff
FAMILIES <- c("gut_cell_atlas", "egozi", "hickey", "gse214852_hio")
CLASSES <- c("Primary-labeled EC", "Other terminal EEC")
sw <- rbind(rd(a[["sweep"]]), rd(a[["hio_sweep"]]))
rf <- rbind(rd(a[["refits"]]), rd(a[["hio_refits"]]))
stopifnot(setequal(unique(sw$family), FAMILIES), setequal(unique(rf$family), FAMILIES))

curve <- sw[sw$class %in% CLASSES & sw$family %in% FAMILIES, ]
on_grid <- abs(curve$cutoff * 50 - round(curve$cutoff * 50)) < 1e-9 | abs(curve$cutoff - CUT) < 1e-12
curve <- curve[on_grid, ]                                       # the 0.02 grid and the cutoff used
refit <- rf[rf$class == "Primary-labeled EC" & rf$family %in% FAMILIES, ]
stopifnot(nrow(refit) == 20L, all(table(refit$family) == 5L), length(unique(refit$equiv_cutoff)) == 5L,
          sum(abs(curve$cutoff - CUT) < 1e-12) == 2L * length(FAMILIES))
out <- rbind(data.frame(kind = "curve", curve[c("family", "class", "cutoff", "post_rate")]),
             data.frame(kind = "refit", family = refit$family, class = refit$class, cutoff = refit$equiv_cutoff,
                        post_rate = refit$post_rate))
write.table(out, a[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
