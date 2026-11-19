# Figure S3C (and Figures 1E, 4B through the shared results): the four-marker score and classification of every
# mature endocrine cell of the six SC-islet studies, with the fitted constants of marker_fit.R.
#
# What it computes
#   Log1p counts per million of the four markers (raw UMI counts over each cell's full-object UMI total), the
#   four-marker score with the fitted standardization, candidate = score >= cutoff and >= 2 of 4 markers detected,
#   and the neighborhood check of each study's 25-nearest-neighbor graph (BH q <= 0.05).
#
# Inputs : constants (fitted_constants.json), metadata (one row per SC endocrine cell), graphs (each study's
#          neighborhood graphs), counts (raw counts of the marker genes), cells (the drawn cells' order and class)
# Outputs: scores (study, class and score, 6 decimals, of every drawn cell), cutoff (the cutoff),
#          classified (per study: SC-SICs, SC-beta cells and other SC endocrine cells, candidates and classified as EC)
#   Rscript score_sc_islet.R constants=... metadata=... graphs=... counts=... cells=... scores=... cutoff=... classified=...
suppressPackageStartupMessages({ library(Matrix); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
STUDY <- c(GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022", GSE185036 = "Du 2022",
           GSE199636 = "Augsornworawat 2023 (nuclei)", GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024")
num <- function(x) { s <- sprintf("%.15g", x); bad <- !is.na(x) & suppressWarnings(as.numeric(s)) != x
  s[bad] <- sprintf("%.17g", x[bad]); s[is.na(x)] <- "NA"; s }
w <- function(df, path) { con <- if (grepl("[.]gz$", path)) gzfile(path, "wb") else file(path, "wb")
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n"); close(con) }

# the neighborhood check (as ec_graph_veto in ec_classifier_lib.R)
ec_graph_veto <- function(candidate, graphs, veto_q = 0.05) {
  local_q <- rep(NA_real_, length(candidate))
  for (g in graphs) {
    if (is.null(g)) next
    is_cand_a <- candidate[g$cols]
    S <- sum(is_cand_a)
    if (S < 2L) next
    hits <- rowSums(matrix(is_cand_a[g$nn], nrow = g$n_a, ncol = g$k))
    foc <- which(is_cand_a)
    p <- phyper(hits[foc] - 1L, S - 1L, (g$n_a - 1L) - (S - 1L), g$k, lower.tail = FALSE)
    local_q[g$cols[foc]] <- p.adjust(p, method = "BH")
  }
  candidate & !is.na(local_q) & local_q <= veto_q
}

# ---- SC cells, in the order of each study's graphs
meta <- read.delim(a[["metadata"]], stringsAsFactors = FALSE, colClasses = "character")
meta$total_umis <- as.numeric(meta$total_umis)
meta$file_order <- seq_len(nrow(meta))
sg <- readRDS(a[["graphs"]])
local({ sg0 <- sg$classifier_gate; for (f in names(sg0)) { idx <- which(meta$family_id == f)
  stopifnot(all(diff(idx) == 1L)); o <- match(sg0[[f]]$cells, meta$project_cell_id[idx]); stopifnot(!anyNA(o), length(o) == length(idx))
  meta[idx, ] <<- meta[idx[o], ] } })
mc <- fromJSON(a[["constants"]])
CORE4 <- mc$ec_core_genes; stopifnot(identical(sort(CORE4), c("FEV", "LMX1A", "SLC18A1", "TPH1")))
sv <- read.delim(gzfile(a[["counts"]]), stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric"))
sv <- sv[sv$gene %in% CORE4, ]
sgenes <- sort(unique(sv$gene), method = "radix"); stopifnot(setequal(sgenes, CORE4))
j <- match(sv$project_cell_id, meta$project_cell_id); stopifnot(!anyNA(j))
expr <- sparseMatrix(i = match(sv$gene, sgenes), j = j, x = log1p(sv$count / meta$total_umis[j] * 1e6),
                     dims = c(length(sgenes), nrow(meta)), dimnames = list(sgenes, meta$project_cell_id))
stopifnot(nrow(meta) == 55120L, all(is.finite(expr@x)))
STUDIES <- sort(unique(meta$family_id)); stopifnot(identical(STUDIES, names(STUDY)))
mgraphs <- sg$marker_rule_veto
stopifnot(identical(names(mgraphs), STUDIES))
for (f in STUDIES) stopifnot(identical(mgraphs[[f]]$cols, which(meta$family_id == f)))

# ---- four-marker score, candidate, neighborhood check
st <- mc$standardisation; rownames(st) <- st$gene
Z <- vapply(CORE4, function(g) (as.numeric(expr[g, ]) - st[g, "mean"]) / st[g, "sd"], numeric(ncol(expr)))
mk_score <- rowMeans(Z)
mk_detected <- as.numeric(Matrix::colSums(expr[CORE4, , drop = FALSE] > 0))
mk_candidate <- mk_score >= mc$ec_core_cutoff & mk_detected >= mc$min_markers_detected
mk_high <- ec_graph_veto(mk_candidate, mgraphs, mc$neighbourhood$q_threshold)

# ---- per-cell table back in metadata (file) order
o <- order(meta$file_order)
cells <- data.frame(project_cell_id = meta$project_cell_id[o], study = unname(STUDY[meta$family_id[o]]), class = meta$target[o],
                    score6 = sprintf("%.6f", mk_score[o]), candidate = mk_candidate[o], classified_ec = mk_high[o],
                    stringsAsFactors = FALSE)
cells$score6[cells$score6 == "-0.000000"] <- "0.000000"
# scores with 6 decimals, rows in the order of the cells table (its cell_id is project_cell_id without the first
# "|"-separated field); that table pools Balboa 2022's SC_EPSILON and SC_GAMMA cells as SC_EPSILON_GAMMA (16 cells),
# and its class is the one written
f3 <- read.delim(gzfile(a[["cells"]]), stringsAsFactors = FALSE, colClasses = "character")
f3 <- f3[f3$dataset %in% names(STUDY), ]                                    # the six studies
suffix <- sub("^[^|]*\\|", "", cells$project_cell_id); stopifnot(!anyDuplicated(suffix), !anyDuplicated(f3$cell_id))
fo <- match(f3$cell_id, suffix); stopifnot(!anyNA(fo), length(fo) == nrow(cells))
dif <- f3$target != cells$class[fo]
stopifnot(all(f3$target[dif] == "SC_EPSILON_GAMMA"), all(cells$class[fo][dif] %in% c("SC_EPSILON", "SC_GAMMA")))
w(data.frame(study = cells$study[fo], class = f3$target, score = num(as.numeric(cells$score6[fo]))), a[["scores"]])
con <- file(a[["cutoff"]], "wb"); writeLines(c("cutoff", num(mc$ec_core_cutoff)), con); close(con)

# ---- classified as EC (score, >= 2 markers and the neighborhood check), per study and class
grp <- list(SC_SIC = function(cl) cl == "SC_SIC", SC_BETA = function(cl) cl == "SC_BETA",
            OTHER_MATURE_SC_ENDOCRINE = function(cl) cl != "SC_SIC")
calls <- do.call(rbind, lapply(names(STUDY), function(f) do.call(rbind, lapply(names(grp), function(g) {
  k <- cells$study == STUDY[[f]] & grp[[g]](cells$class)
  data.frame(held_out_study = f, study = STUDY[[f]], class = g, n = sum(k), candidates = sum(cells$candidate[k]),
             classified_ec = sum(cells$classified_ec[k]), pct = round(100 * mean(cells$classified_ec[k]), 2),
             pct_candidates = round(100 * mean(cells$candidate[k]), 2))
}))))
w(calls, a[["classified"]])
cat(sprintf("cutoff %.17g; %d cells scored\n", mc$ec_core_cutoff, nrow(cells)))
