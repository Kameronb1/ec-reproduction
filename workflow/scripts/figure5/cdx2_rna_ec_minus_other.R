# CDX2 RNA, EC minus other endocrine cells, in the five studies of Figure S17E (its RNA axis: the same contrast as its
# chromatin axis). Same pseudobulk profiles and thinning as Figure 5A's regulator RNA (raw-count sums of cells with primary
# study labels, one profile per sample and cell type; profiles above 569,504 reads binomially thinned to that depth, 20
# draws, CDX2 log2(1 + CPM) averaged over draws). Thinning only CDX2 and the rest of the library separately is the same
# draw as thinning every gene, so the CPM denominator stays exact. Profiles that Figure 5A already drew keep 5A's values
# (tf_expression_by_sample.tsv), so the EC levels equal 5A's; the populations 5A does not draw (SC-delta, Hickey 2023's and
# McCauley 2023's other enteroendocrine cells, McCauley 2023's EC) are drawn here.
# Contrast: per sample, EC minus each other endocrine population with at least 20 cells (EC also at least 20), averaged over
# those populations; the study value is the mean over its samples. Other endocrine populations, as in each study's chromatin
# contrast: SC-islet studies SC-beta, SC-alpha, SC-delta; Singh 2024 D, X and T4 cells; Hickey 2023 every other
# enteroendocrine cell; McCauley 2023 every other enteroendocrine cell of the Endocrine_1 subset (no chromatin data).
#
#   Rscript cdx2_rna_ec_minus_other.R sc_sums= sc_groups= gut_sums= gut_groups= induced_sums= induced_groups= hio_sums=
#           hio_groups= by_sample=<tf_expression_by_sample.tsv> out=<per study> by_profile=<per profile>

a <- local({ x <- commandArgs(TRUE); setNames(sub("^[^=]+=", "", x), sub("=.*$", "", x)) })
DEPTH <- 569504; REPS <- 20L; SEED <- 20260925L

STUDY <- list(
  gse199636 = list(src = "sc", family = "GSE199636", ec = "SC_SIC", other = c("SC_BETA", "SC_ALPHA", "SC_DELTA")),
  zhu = list(src = "sc", family = "GSE202497", ec = "SC_SIC", other = c("SC_BETA", "SC_ALPHA", "SC_DELTA")),
  singh = list(src = "induced", family = "SF_INTESTINAL_GSE238274", ec = "EC", other = c("D", "X", "T4")),
  hickey = list(src = "gut", family = "SF_HICKEY2023_ADULT_EEC", ec = "POSITIVE", other = "NEGATIVE"),
  mccauley = list(src = "hio", family = "SF_GSE214852_HIO_GRAFT", ec = "EC", other = "OTHER_EEC"))

# the header and the CDX2 row of a gene x profile sums table
cdx2_row <- function(f) {
  con <- gzfile(f, "r"); on.exit(close(con))
  hdr <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]]
  repeat { l <- readLines(con, 5000); if (!length(l)) stop("no CDX2 row in ", f)
    k <- which(startsWith(l, "CDX2\t")); if (length(k)) break }
  setNames(as.numeric(strsplit(l[k[1]], "\t", fixed = TRUE)[[1]][-1]), hdr[-1])
}
P <- do.call(rbind, lapply(unique(vapply(STUDY, `[[`, "", "src")), function(s) {
  g <- read.delim(a[[paste0(s, "_groups")]], stringsAsFactors = FALSE, check.names = FALSE)
  x <- cdx2_row(a[[paste0(s, "_sums")]])
  g$cdx2 <- unname(x[g$group]); stopifnot(!anyNA(g$cdx2)); g }))
P <- P[P$family_id %in% vapply(STUDY, `[[`, "", "family"), ]
set.seed(SEED)
P$level <- vapply(seq_len(nrow(P)), function(j) {
  c0 <- P$cdx2[j]; lib <- P$lib[j]
  if (lib < DEPTH) return(log2(1 + c0 / lib * 1e6))
  mean(replicate(REPS, { k <- rbinom(1, c0, DEPTH / lib); r <- rbinom(1, lib - c0, DEPTH / lib); log2(1 + k / (k + r) * 1e6) }))
}, 0)
TR <- read.delim(a[["by_sample"]], stringsAsFactors = FALSE, check.names = FALSE); TR <- TR[TR$tf == "CDX2", ]
m <- match(paste(P$family_id, P$unit, P$class), paste(TR$family, TR$unit, TR$population))
P$source <- ifelse(is.na(m), "drawn here", "Figure 5A table"); P$level[!is.na(m)] <- TR$mean_log2cpm[m[!is.na(m)]]
write.table(P[, c("group", "family_id", "unit", "class", "n_cells", "lib", "cdx2", "level", "source")], a[["by_profile"]],
            sep = "\t", quote = FALSE, row.names = FALSE)

OUT <- do.call(rbind, lapply(names(STUDY), function(k) { s <- STUDY[[k]]; d <- P[P$family_id == s$family & P$n_cells >= 20, ]
  per <- do.call(rbind, lapply(unique(d$unit[d$class == s$ec]), function(u) { e <- d[d$unit == u & d$class == s$ec, ]; o <- d[d$unit == u & d$class %in% s$other, ]
    if (!nrow(o)) return(NULL)
    data.frame(unit = u, ec = e$level, other = mean(o$level), diff = mean(e$level - o$level), others = paste(o$class, collapse = "+")) }))
  data.frame(cohort = k, samples = nrow(per), ec_rna = mean(per$ec), other_rna = mean(per$other), ec_minus_other = mean(per$diff),
             min_sample = min(per$diff), max_sample = max(per$diff), others_used = paste(unique(per$others), collapse = "; ")) }))
write.table(OUT, a[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
