# Figure S5A: how the 20-gene classifier classifies a cell, worked through for the Elmentaite 2021 model
# (trained on Egozi 2021 and Hickey 2023). The model is refitted with the classifier's functions
# (classifier_lib/) and must equal the weights, training effects and cutoff written by classifier.R
# (1e-12), and Fig. 2D's per-study values (1e-9).
#
#   genes    per gene: Egozi 2021 and Hickey 2023 differences (EC minus other terminal EEC, donor means),
#            their median, the weight (median / sum of medians), and the z-score mean and SD
#   mixture  the score histogram of cells with >= 5 of 20 genes detected (studies weighted equally,
#            bins of 0.1), the two fitted Gaussian components and their mixture (600 points), and the
#            mixture parameters (one row per parameter)
#   hist     the drawn score histogram (bin start, bin end, share of cells), 17 significant digits
#   mix      the drawn mixture: component weights and means, shared SD and cutoff (one row), 17 significant digits
#
# Inputs : analysis (training corpus, RDS), weights, predictions (classifier.R), d_genes (Fig. 2D gene table)
# Outputs: genes, mixture, hist, mix
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
pkg <- file.path(HERE, "classifier_lib")
for (f in c("figure2_postvalidation_core.R", "figure2_postvalidation_analysis.R"))
  source(file.path(pkg, "vendor/base/R", f))
source(file.path(pkg, "R/candidate_panel_core.R")); source(file.path(pkg, "R/population_policy.R"))
# Matrix is attached before any fit (see cp_fit_standardizer)
suppressPackageStartupMessages(library(Matrix)); stopifnot("package:Matrix" %in% search())
x <- readRDS(a[["analysis"]]); classifier_assert_training_data(x)
m <- x$metadata; labels <- classifier_preflight_labels(m)
contract <- cp_read_json(file.path(pkg, "config/candidate_panel_contract.json"))
g <- contract$gmm
gmm <- list(deterministic_initial_quantile_pairs = g$initial_quantile_pairs, max_iterations = as.integer(g$max_iterations),
            tolerance = as.numeric(g$tolerance), variance_floor = as.numeric(g$variance_floor),
            minimum_component_weight = as.numeric(g$minimum_component_weight), minimum_component_sd = as.numeric(g$minimum_component_sd),
            minimum_mean_separation_pooled_sd = as.numeric(g$minimum_mean_separation_pooled_sd))
ELM <- "SF_ELMENTAITE2021_GUT_ATLAS"; EGO <- "SF_EGOZI2021_DEVELOPMENTAL_EEC"; HIC <- "SF_HICKEY2023_ADULT_EEC"
TRAIN_FAM <- c(EGO, HIC)
w <- read.delim(a[["weights"]], stringsAsFactors = FALSE); w <- w[w$study == ELM & w$panel_size == 20, ]
stopifnot(nrow(w) == 20L, identical(w$rank, 1:20), identical(w$trained_on, rep(paste(TRAIN_FAM, collapse = ";"), 20)))
D <- read.delim(a[["d_genes"]], stringsAsFactors = FALSE)
GENES <- unique(D$gene_symbol[order(D$panel_rank)]); stopifnot(identical(GENES, w$gene_symbol))
pred <- read.delim(gzfile(a[["predictions"]]), stringsAsFactors = FALSE)
THR <- unique(pred$threshold[pred$study == ELM & pred$panel_size == 20]); stopifnot(length(THR) == 1L)
REQ <- unique(pred$detection_requirement[pred$study == ELM & pred$panel_size == 20]); stopifnot(identical(as.integer(REQ), 5L))

# refit the Elmentaite 2021 size-20 model; it must reproduce the stored model
train <- m$family_id %in% TRAIN_FAM
fit <- cp_fit_panel_model(x$normalized, x$genes, m, labels, train, GENES, gmm, detection_fraction = 0.25)
mod <- fit$model
stopifnot(max(abs(mod$gene_weights - w$weight)) < 1e-12, max(abs(mod$training_family_effects - w$training_effect)) < 1e-12,
          identical(mod$gmm$status, "PASS"), abs(mod$gmm$threshold - THR) < 1e-12, mod$detection_requirement == REQ,
          identical(mod$training_families, sort(TRAIN_FAM)))
# z-scores: mean and SD over all training-study cells, each study weighted equally
bw <- cp_family_balanced_weights(m, train)
stopifnot(all(abs(tapply(bw, m$family_id[train], sum) - .5) < 1e-12))
xm <- as.matrix(x$normalized[match(GENES, rownames(x$normalized)), train])
mu <- as.numeric(xm %*% bw); sdv <- sqrt(pmax(as.numeric((xm^2) %*% bw) - mu^2, 0))
stopifnot(max(abs(mu - mod$standardizer$training_mean)) < 1e-10, max(abs(sdv - mod$standardizer$training_sd)) < 1e-10)
# weights: per-study donor-mean differences, median, normalized to sum to 1
u <- cp_unit_class_means(x$normalized, m, labels, train, GENES, 1L)
eff <- sapply(TRAIN_FAM, function(f) { b <- u$metadata$family_id == f
  rowMeans(u$values[, b & u$metadata$label == "POSITIVE", drop = FALSE]) - rowMeans(u$values[, b & u$metadata$label == "NEGATIVE", drop = FALSE]) })
dD <- sapply(TRAIN_FAM, function(f) D$effect[D$family_id == f][match(GENES, D$gene_symbol[D$family_id == f])])
stopifnot(max(abs(eff - dD)) < 1e-9, all(D$label_source[D$family_id %in% TRAIN_FAM] == "AUTHOR"),
          max(abs(apply(dD, 1, median) - w$training_effect)) < 1e-9, all(w$training_effect > 0),
          max(abs(w$training_effect / sum(w$training_effect) - w$weight)) < 1e-12)
# score = weighted mean of the 20 z-scores
z <- (as.matrix(x$normalized[match(GENES, rownames(x$normalized)), ]) - mu) / sdv
stopifnot(max(abs(as.numeric(crossprod(w$weight, z)) - fit$score)) < 1e-9, abs(sum(w$weight) - 1) < 1e-12)
# cutoff: where the two weighted components of the equal-variance mixture cross
gc <- train & fit$detected_genes >= REQ & is.finite(fit$score)
gw <- cp_family_balanced_weights(m, gc)
G <- mod$gmm
cross <- function(v) G$lower_weight * dnorm(v, G$lower_mean, G$shared_sd) - G$upper_weight * dnorm(v, G$upper_mean, G$shared_sd)
stopifnot(abs(cross(G$threshold)) < 1e-10, G$threshold > G$lower_mean, G$threshold < G$upper_mean)
stopifnot(identical(c(sprintf("%.2f", c(dD[1, EGO], dD[1, HIC], w$training_effect[1])), sprintf("%.1f", sum(w$training_effect)),
                      sprintf("%.1f%%", 100 * w$weight[1]), sprintf("%.2f", THR)), c("4.73", "4.10", "4.41", "44.9", "9.8%", "0.50")))

write.table(data.frame(gene = GENES, egozi_2021 = dD[, EGO], hickey_2023 = dD[, HIC], median = w$training_effect, weight = w$weight,
                       z_mean = mu, z_sd = sdv), a[["genes"]], sep = "\t", quote = TRUE, row.names = FALSE)
BW <- .1; sc <- fit$score[gc]
brk <- seq(floor(min(sc) / BW) * BW, ceiling(max(sc) / BW) * BW + 1e-9, BW)
h <- as.numeric(tapply(gw, cut(sc, brk, right = FALSE, include.lowest = TRUE), sum)); h[is.na(h)] <- 0
stopifnot(abs(sum(h) - 1) < 1e-9)
xs <- seq(min(brk), max(brk), length.out = 600)
lo <- BW * G$lower_weight * dnorm(xs, G$lower_mean, G$shared_sd); up <- BW * G$upper_weight * dnorm(xs, G$upper_mean, G$shared_sd)
out <- rbind(data.frame(part = "histogram", x = head(brk, -1), x_end = brk[-1], value = h),
             data.frame(part = "lower_component", x = xs, x_end = NA, value = lo),
             data.frame(part = "upper_component", x = xs, x_end = NA, value = up),
             data.frame(part = "mixture_fit", x = xs, x_end = NA, value = lo + up),
             data.frame(part = c("lower_mean", "upper_mean", "lower_weight", "upper_weight", "shared_sd", "cutoff", "cells", "cells_egozi", "cells_hickey"),
                        x = NA, x_end = NA, value = c(G$lower_mean, G$upper_mean, G$lower_weight, G$upper_weight, G$shared_sd, G$threshold,
                                                      sum(gc), sum(gc & m$family_id == EGO), sum(gc & m$family_id == HIC))))
write.table(out, a[["mixture"]], sep = "\t", quote = TRUE, row.names = FALSE)

# the tables S5A (iii) draws: histogram bins and the six mixture parameters, 17 significant digits
put <- function(d, name) write.table(data.frame(lapply(d, sprintf, fmt = "%.17g")), a[[name]], sep = "\t", quote = FALSE,
                                     row.names = FALSE)
put(data.frame(lo = head(brk, -1), hi = brk[-1], h = h), "hist")
put(data.frame(lower_weight = G$lower_weight, lower_mean = G$lower_mean, upper_weight = G$upper_weight,
               upper_mean = G$upper_mean, shared_sd = G$shared_sd, threshold = G$threshold), "mix")
