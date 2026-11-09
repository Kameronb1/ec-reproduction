# Compare a drawn panel with the published one and write a JSON report.
#
#   Rscript compare_image.R <observed.png> <expected.png> <report.json> <exact|tolerance>
#
# EXACT             every pixel is identical. Expected on the reference setup: R 4.6.1 on Windows, cairo PNG
#                   device, Arial, the package versions in envs/r_packages.tsv.
# WITHIN_TOLERANCE  not identical, but the same picture: after averaging both images over 24 x 24 pixel blocks
#                   (1/25 in), no block differs by more than 64 of 255 levels in any color channel. Other systems
#                   draw text and line edges with different anti-aliasing, font hinting and dash patterns; the same
#                   panels drawn with another renderer (Windows GDI instead of cairo) differ by up to 49 levels.
#                   A moved, missing or recolored mark changes its blocks by more: a point moved by 15 percentage
#                   points differs by 188 levels, a cutoff line moved by 0.1 by 99. A panel that places marks by
#                   measuring its text can shift as a whole where fonts measure differently, and then fails.
# DIFFERENT         anything else, including a different image size.
#
# With "exact", only EXACT passes; with "tolerance", EXACT and WITHIN_TOLERANCE pass. A failing comparison exits
# non-zero, so it stops the workflow.

args <- commandArgs(TRUE)
stopifnot(length(args) == 4, args[4] %in% c("exact", "tolerance"))
BLOCK <- 24; LEVELS <- 64

rgb_of <- function(path) {
  a <- png::readPNG(path)
  if (length(dim(a)) == 2) a <- array(a, c(dim(a), 3))
  if (dim(a)[3] == 2) a <- a[, , c(1, 1, 1), drop = FALSE]
  if (dim(a)[3] == 4) {                                          # composite any transparency on white
    alpha <- a[, , 4]
    a <- a[, , 1:3] * as.vector(alpha) + (1 - as.vector(alpha))
  }
  round(a[, , 1:3] * 255)
}
blocks <- function(m) {                                          # mean of each BLOCK x BLOCK block of a matrix
  gr <- (seq_len(nrow(m)) - 1) %/% BLOCK; gc <- (seq_len(ncol(m)) - 1) %/% BLOCK
  rowsum(t(rowsum(m, gr)), gc) / outer(tabulate(gc + 1), tabulate(gr + 1))
}

obs <- rgb_of(args[1]); exp <- rgb_of(args[2])
report <- list(observed = args[1], expected = args[2], width = dim(exp)[2], height = dim(exp)[1])
if (!identical(dim(obs), dim(exp))) {
  report$status <- "DIFFERENT"
  report$detail <- sprintf("size %d x %d, expected %d x %d", dim(obs)[2], dim(obs)[1], dim(exp)[2], dim(exp)[1])
} else {
  d <- abs(obs - exp)
  pix <- apply(d, c(1, 2), max)
  bd <- do.call(pmax, lapply(1:3, function(k) abs(blocks(obs[, , k]) - blocks(exp[, , k]))))
  report$pixels_differing <- sum(pix > 0)
  report$pixels_differing_share <- mean(pix > 0)
  report$max_pixel_difference <- max(pix)
  report$max_block_difference <- max(bd)
  report$status <- if (max(pix) == 0) "EXACT" else if (max(bd) <= LEVELS) "WITHIN_TOLERANCE" else "DIFFERENT"
}
report$required <- args[4]
passed <- report$status == "EXACT" || (args[4] == "tolerance" && report$status == "WITHIN_TOLERANCE")
dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(report, args[3], auto_unbox = TRUE, pretty = TRUE, digits = NA)
# Snakemake deletes the report of a failed job, so its numbers are printed here too.
shown <- report[setdiff(names(report), c("observed", "expected", "status", "required"))]
if (!passed) stop(args[1], " differs from ", args[2], ": ", report$status, ", ",
                  paste(names(shown), unlist(shown), sep = " = ", collapse = ", "), call. = FALSE)
