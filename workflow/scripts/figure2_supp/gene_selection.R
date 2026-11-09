# Figure S4A-B: how stable the 20 genes are when each gene-ranking dataset is left out, and whether
# that recurrence is more than chance. Both come from deposited tables of Figure 2's leave-one-study-out
# run: its seven leave-one-study-out gene rankings and their recurrence null (128 exact study-block
# label swaps, 200 within-donor label swaps, each re-ranking all seven leave-one-study-out lists).
#
#   S4A  for each gene (the 20 of Fig. 2D, then the others reaching any top-20 list) and each left-out
#          study: in that study's top 20, at ranks 21-30, or below 30
#   S4B  genes in all seven top-20 lists for every label-swap replicate (size 20, completed replicates);
#          with the real labels (observed: genes in all seven lists of size 10, 20 and 30, and how many
#          distinct genes the seven lists hold); and per swap level and list size the null mean, 97.5th
#          percentile and one-sided p of the real count (summary). The study-block null enumerates all
#          2^7 = 128 assignments, so its p is exact: replicates with at least the real count / 128. The
#          within-donor null is a Monte Carlo sample: (1 + replicates with at least the real count) / (200 + 1).
#
# Inputs : panels (nested leave-one-study-out panels), d_genes (Fig. 2D gene table), replicates (recurrence null)
# Outputs: grid, null, observed, summary
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))

pan <- read.delim(gzfile(a[["panels"]]), stringsAsFactors = FALSE, check.names = FALSE)
p20 <- pan[pan$panel_id == "NESTED_SIZE_20", ]; p30 <- pan[pan$panel_id == "NESTED_SIZE_30", ]
dt <- read.delim(a[["d_genes"]], stringsAsFactors = FALSE)
D <- unique(dt$gene_symbol[order(dt$panel_rank)]); stopifnot(length(D) == 20L)
genes <- unique(c(D, p20$gene_symbol))
cnt <- sapply(genes, function(g) sum(p20$gene_symbol == g))
# genes in order: the 20 of Fig. 2D by number of lists then rank, then the others by number of lists then name
ord <- c(D[order(-cnt[D], match(D, D))], setdiff(genes, D)[order(-cnt[setdiff(genes, D)], setdiff(genes, D))])
# left-out studies in Fig. 2D's order: studies labeled by the four-marker classifier, by display name, then the
# studies with primary study labels
FOLD_ORD <- c("SF_GSE185224_ADULT_REGIONAL_INTESTINE", "SF_ELMENTAITE2020_DEVGUT", "SF_GSE125970_ADULT_INTESTINE",
              "SF_YU2021_ENDODERM", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC", "SF_EGOZI2021_DEVELOPMENTAL_EEC")
stopifnot(setequal(FOLD_ORD, unique(p20$outer_family)))
grid_df <- expand.grid(gene = ord, fold = FOLD_ORD, stringsAsFactors = FALSE)
key <- paste(grid_df$gene, grid_df$fold)
grid_df$state <- ifelse(key %in% paste(p20$gene_symbol, p20$outer_family), "Top 20",
                 ifelse(key %in% paste(p30$gene_symbol, p30$outer_family), "Ranks 21-30", "Below 30"))
write.table(grid_df, a[["grid"]], sep = "\t", quote = TRUE, row.names = FALSE)   # the drawn key prints 21-30 with an en dash

rr <- read.delim(a[["replicates"]], stringsAsFactors = FALSE)
rep <- rr[rr$panel_size == 20 & rr$status == "PASS", ]
stopifnot(sum(rep$swap_level == "FAMILY_BLOCK") == 128L, sum(rep$swap_level == "UNIT") == 200L)
write.table(rep[c("swap_level", "replicate", "genes_in_every_fold")], a[["null"]], sep = "\t", quote = TRUE, row.names = FALSE)

# real labels: genes in all seven leave-one-study-out lists of each size (names in byte order, whatever the locale)
SIZES <- c(10L, 20L, 30L)
observed <- do.call(rbind, lapply(SIZES, function(k) {
  p <- pan[pan$panel_id == paste0("NESTED_SIZE_", k), ]
  sets <- split(p$gene_symbol, p$outer_family); stopifnot(length(sets) == 7L, all(lengths(sets) == k))
  shared <- sort(Reduce(intersect, sets), method = "radix")
  data.frame(panel_size = k, genes_in_every_fold = length(shared), distinct_genes = length(unique(unlist(sets))),
             every_fold_genes = paste(shared, collapse = "; "))
}))
write.table(observed, a[["observed"]], sep = "\t", quote = FALSE, row.names = FALSE)

# null summary per swap level and list size, over every attempted replicate (failed ones count against the p bounds)
stopifnot(all(table(rr$swap_level, rr$panel_size)["FAMILY_BLOCK", ] == 2^7))   # all study-block assignments: exact null
summary <- do.call(rbind, lapply(split(rr, paste(rr$swap_level, rr$panel_size)), function(z) {
  good <- z$status == "PASS"; v <- z$genes_in_every_fold[good]; bad <- sum(!good)
  obs <- observed$genes_in_every_fold[match(z$panel_size[1], observed$panel_size)]; ex <- sum(v >= obs)
  unit <- z$swap_level[1] == "UNIT"; pseudo <- as.integer(unit); den <- nrow(z) + pseudo   # plus-one for Monte Carlo only
  data.frame(swap_level = z$swap_level[1], panel_size = z$panel_size[1], observed = obs, attempted = nrow(z),
             successful = sum(good), failed = bad, null_mean = if (length(v)) mean(v) else NA_real_,
             null_q975 = if (length(v)) unname(quantile(v, .975)) else NA_real_, exceedances = ex,
             empirical_p = if (bad == 0L) (pseudo + ex) / den else NA_real_, p_lower = (pseudo + ex) / den,
             p_upper = (pseudo + ex + bad) / den,
             tail_method = if (unit) "MONTE_CARLO_PLUS_ONE" else "EXACT_ALL_FAMILY_ASSIGNMENTS")
}))
summary <- summary[order(summary$swap_level, summary$panel_size, method = "radix"), ]
write.table(summary, a[["summary"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
