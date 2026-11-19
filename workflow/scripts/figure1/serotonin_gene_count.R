# Figure 1D (last column) and Supplementary Figure S1A: how many of the three serotonin genes each cell carries, and the
# islet hormone transcripts in each cell group.
# Cells, groups, Crohn's exclusion and detection rule are those of serotonin_codetection.R (Figure 1D): detection is
# >= 0.5 SoupX-corrected count (unrounded).
#   out_count   per study and cell group: share of cells with 0, 1, 2 or 3 of TPH1, DDC and SLC18A1 detected, and the
#               TPH1 share; drawn is taken from the Figure 1D table (groups with >= 30 cells). The figure draws, per
#               cell group, the mean of the drawn studies' shares.
#   out_groups  per study and cell group: the share of cells with INS, GCG and SST detected (each >= 0.5 SoupX-corrected
#               count) and with a hormone of another islet cell type (beta cells: GCG or SST; alpha cells: INS or SST;
#               delta cells: INS or GCG; intestinal EC, HIO graft EC and SC-SIC: any of the three), in all cells of the
#               group (base "all") and in its cells with TPH1 detected (base "tph1").
# Groups: intestinal EC (primary intestine EC label; table group "Gut EC"), HIO graft EC (McCauley 2023 Endocrine_1),
# intestinal other EEC (primary intestine "Other EEC" label; table group "Gut, other EEC"), SC-SIC, SC-beta, SC-alpha,
# SC-delta, adult beta, alpha, delta.
# The 7 Crohn's disease donors are left out; a cell without corrected counts is left out.
# Inputs: cells        per-study display cells (study, context, sample, published label, class, uncorrected TPH1/DDC/SLC18A1);
#                      the context "Native intestine" is the primary intestine
#         sc           SoupX raw and corrected counts, SC-islet cells
#         gut_ec       SoupX raw and corrected counts, intestinal EC cells (column corrected = SoupX applied)
#         gut_other    SoupX raw and corrected counts, primary intestine other EEC cells
#         adult        SoupX raw and corrected counts, adult islet endocrine cells
#         codetection  the Figure 1D table written by serotonin_codetection.R (group "HIO EC" is read as "HIO graft EC");
#                      supplies drawn and the checks below
#   Rscript serotonin_gene_count.R cells=... sc=... gut_ec=... gut_other=... adult=... codetection=<serotonin_codetection.tsv>
#           out_count=<serotonin_gene_count.tsv> out_groups=<hormones_by_group.tsv>
args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
stopifnot(all(c("cells", "sc", "gut_ec", "gut_other", "adult", "codetection", "out_count", "out_groups") %in% names(arg)))

KIT <- c("TPH1", "DDC", "SLC18A1"); HORM <- c("INS", "GCG", "SST")
CROHNS <- c("T017", "T019", "T176", "T189", "T197", "T202", "T203")
check <- function(name, cond) if (!isTRUE(cond)) stop("check failed: ", name, call. = FALSE)

# ---- cells and groups (as Figure 1D) ---------------------------------------------------------------------------------
cl <- read.delim(gzfile(arg[["cells"]]), check.names = FALSE, stringsAsFactors = FALSE)
cl$study[cl$study == "Elmentaite 2020 + Elmentaite 2021"] <- "Elmentaite 2020 + 2021"
cl <- cl[!cl$sample %in% CROHNS, ]
nat <- cl$context == "Native intestine"
cl$grp <- NA_character_
cl$grp[nat & cl$class == "EC"] <- "Gut EC"
cl$grp[cl$study == "McCauley 2023" & cl$author_label == "Endocrine_1"] <- "HIO graft EC"
cl$grp[nat & cl$class == "Other EEC"] <- "Gut, other EEC"
for (k in c("SC-SIC", "SC-beta", "SC-alpha", "SC-delta")) cl$grp[cl$context == "SC-islet" & cl$class == k] <- k
for (k in c("Beta", "Alpha", "Delta")) cl$grp[cl$context == "Adult islet" & cl$class == k] <- paste("Adult", tolower(k))
cl <- cl[!is.na(cl$grp), ]

# ---- corrected counts and detection ----------------------------------------------------------------------------------
rd <- function(k) read.delim(gzfile(arg[[k]]), stringsAsFactors = FALSE, check.names = FALSE)
sx <- rd("sc"); ge <- rd("gut_ec"); go <- rd("gut_other"); ad <- rd("adult"); ge <- ge[ge$corrected, ]; go <- go[go$corrected, ]
cols <- paste0(c(KIT, HORM), "_soupx")
CORR <- rbind(sx[, c("cell", cols)], ge[, c("cell", cols)], go[, c("cell", cols)], ad[, c("cell", cols)])
check("corrected cells unique", !anyDuplicated(CORR$cell))
cl <- cl[cl$cell_id %in% CORR$cell, ]; m <- match(cl$cell_id, CORR$cell)
for (g in c(KIT, HORM)) cl[[paste0(g, "_on")]] <- CORR[[paste0(g, "_soupx")]][m] >= 0.5
cl$n_kit <- cl$TPH1_on + cl$DDC_on + cl$SLC18A1_on

# ---- Figure 1D last column: 0, 1, 2 or 3 of the genes ----------------------------------------------------------------
PD <- read.delim(arg[["codetection"]], stringsAsFactors = FALSE); PD$grp[PD$grp == "HIO EC"] <- "HIO graft EC"
D <- do.call(rbind, lapply(split(cl, list(cl$study, cl$grp), drop = TRUE), function(e) data.frame(study = e$study[1], grp = e$grp[1],
  cells = nrow(e), n0 = mean(e$n_kit == 0), n1 = mean(e$n_kit == 1), n2 = mean(e$n_kit == 2), n3 = mean(e$n_kit == 3), TPH1 = mean(e$TPH1_on))))
k <- match(paste(D$study, D$grp), paste(PD$study, PD$grp))
check("every group in the Figure 1D table", !anyNA(k) && nrow(D) == nrow(PD))
check("cells equal the Figure 1D table", all(D$cells == PD$cells[k]))
check("3 of 3 equals the Figure 1D all_three", all(abs(D$n3 - PD$all_three[k]) < 1e-12))
check("TPH1 share equals the Figure 1D table", all(abs(D$TPH1 - PD$TPH1[k]) < 1e-12))
D$drawn <- PD$drawn[k]
write.table(D, arg[["out_count"]], sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Supplementary Figure S1A: islet hormones in every cell group ----------------------------------------------------
OTHER <- list(beta = c("GCG", "SST"), alpha = c("INS", "SST"), delta = c("INS", "GCG"), ec = c("INS", "GCG", "SST"))
kind <- function(g) if (grepl("beta", g, ignore.case = TRUE)) "beta" else if (grepl("alpha", g, ignore.case = TRUE)) "alpha" else if (grepl("delta", g, ignore.case = TRUE)) "delta" else "ec"
cl$other_on <- FALSE
for (g in unique(cl$grp)) { k <- cl$grp == g; o <- OTHER[[kind(g)]]
  cl$other_on[k] <- Reduce(`|`, lapply(o, function(h) cl[[paste0(h, "_on")]][k])) }
G <- do.call(rbind, lapply(c("all", "tph1"), function(base) {
  x <- if (base == "all") cl else cl[cl$TPH1_on, ]
  do.call(rbind, lapply(split(x, list(x$study, x$grp), drop = TRUE), function(e) data.frame(study = e$study[1], grp = e$grp[1], base = base,
    other_hormones = paste(OTHER[[kind(e$grp[1])]], collapse = " or "), cells = nrow(e),
    n_INS = sum(e$INS_on), n_GCG = sum(e$GCG_on), n_SST = sum(e$SST_on), n_other = sum(e$other_on),
    INS = mean(e$INS_on), GCG = mean(e$GCG_on), SST = mean(e$SST_on), other = mean(e$other_on)))) }))
kt <- G$base == "tph1"; ke <- match(paste(G$study[kt], G$grp[kt]), paste(PD$study, PD$grp))
check("TPH1+ cells equal the Figure 1D TPH1 share x cells", !anyNA(ke) && all(G$cells[kt] == round(PD$TPH1[ke] * PD$cells[ke])))
write.table(G, arg[["out_groups"]], sep = "\t", quote = FALSE, row.names = FALSE)
