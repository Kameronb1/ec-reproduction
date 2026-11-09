# Figure S12D: STAC across SC-SIC, EC cells, islet cells and brain neurons, one point per sample.
#
# What it computes. For every sample (donor, organoid or sequenced sample) with >= 10 cells of a group: the % of the
# group's cells with >= 1 STAC count, and STAC ln(1 + CPM) of the sample's pooled counts. Groups: primary intestinal EC
# (primary study labels, Elmentaite 2020 Crohn's donors removed), HIO graft EC, induced EC, SC-SIC, SC-beta, SC-alpha,
# SC-delta, adult islet beta, alpha and delta (beta_major + beta_minor summed per donor; Bandesh 2026 Proliferating_Alpha
# pooled with Alpha per donor, as in Figure 5A), adult brain serotonin and
# dopamine neurons (Siletti 2023, single nuclei) and fetal brain serotonin neurons (Braun 2023, whole cells).
# Two tables: the per-sample table for all samples, and the table the figure draws, in which the nuclei flag is set per
# study (single-nucleus studies: GSE199636, Hickey 2023, Siletti 2023) and the adult islet rows are adults only (the
# donor listed as excluded in islet_roster, HPAP034 aged 13, is removed).
#
# Inputs (name=path): <src>_sums, <src>_detected, <src>_groups for src in gut, hio, ind, sc, islet (per-sample
#   pseudobulk sums, detected-cell counts and group tables); gut_roster (JSON, excluded Crohn's donors);
#   siletti_sums, siletti_detected, siletti_groups, braun_sums, braun_detected, braun_groups (per-donor brain
#   pseudobulk); islet_roster (JSON, adult islet datasets and excluded donors).
# Outputs (name=path): by_sample (all samples), drawn (adult islet rows adults only, nuclei flag per study).
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
GENE <- "STAC"; FLOOR <- 10L

# ---- non-brain samples ----
rd1 <- function(f) { m <- as.data.table(read.delim(gzfile(f), check.names = FALSE)); stopifnot(GENE %in% m$gene); unlist(m[gene == GENE, -1]) }
src <- c("gut", "hio", "ind", "sc", "islet")
U <- rbindlist(lapply(src, function(k) { g <- fread(A[[paste0(k, "_groups")]]); s <- rd1(A[[paste0(k, "_sums")]]); d <- rd1(A[[paste0(k, "_detected")]])
  stopifnot(identical(names(s), g$group), identical(names(d), g$group)); cbind(g, src = k, umi = as.numeric(s), det = as.numeric(d)) }), fill = TRUE)
crohns <- fromJSON(A[["gut_roster"]])$excluded_units
U <- U[!(src == "gut" & family_id == "SF_ELMENTAITE2020_DEVGUT" & unit %in% crohns)]
U[, grp := fcase(src == "gut" & class == "POSITIVE", "Primary gut EC", src == "hio" & class == "EC", "HIO EC", src == "ind" & class == "EC", "Induced EC",
                 src == "sc" & class == "SC_SIC", "SC-SIC", src == "sc" & class == "SC_BETA", "SC-beta", src == "sc" & class == "SC_ALPHA", "SC-alpha",
                 src == "sc" & class == "SC_DELTA", "SC-delta", src == "islet" & class %in% c("Beta", "beta_major", "beta_minor"), "Adult islet beta",
                 src == "islet" & class %in% c("Alpha", "alpha", "Proliferating_Alpha"), "Adult islet alpha", src == "islet" & class %in% c("Delta", "delta"), "Adult islet delta", default = NA_character_)]
U <- U[!is.na(grp)]
U <- U[, .(n_cells = sum(n_cells), lib = sum(lib), umi = sum(umi), det = sum(det)), by = .(grp, study = family_id, unit)]   # beta_major + beta_minor per donor
U[, nuclei := study == "GSE199636"]

# ---- brain samples, per donor (pseudobulk of the cluster's cells) ----
B <- rbindlist(lapply(c("siletti", "braun"), function(k) { g <- fread(A[[paste0(k, "_groups")]]); s <- rd1(A[[paste0(k, "_sums")]]); d <- rd1(A[[paste0(k, "_detected")]])
  g[, .(grp, study, unit, n_cells, lib, umi = as.numeric(s[group]), det = as.integer(d[group]), nuclei)] }))
stopifnot(!anyNA(B$umi), !anyNA(B$det))
A_ <- rbind(U, B)[n_cells >= FLOOR]
A_[, `:=`(pct = 100 * det / n_cells, lcp = log1p(umi / lib * 1e6))]
fwrite(A_, A[["by_sample"]], sep = "\t")

# ---- the drawn table: nuclei flag per study, adult islet rows adults only ----
NUC_STUDY <- c("GSE199636", "SF_HICKEY2023_ADULT_EEC", "Siletti 2023")
CELL_STUDY <- c("GSE151117", "GSE167880", "GSE185036", "GSE202497", "GSE226346", "SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_ELMENTAITE2020_DEVGUT",
                "SF_ELMENTAITE2021_GUT_ATLAS", "SF_GSE214852_HIO_GRAFT", "SF_INTESTINAL_GSE238274", "GSE221156", "GSE148073", "Braun 2023")
drawn_tab <- fread(A[["by_sample"]])
stopifnot(all(drawn_tab$study %in% c(NUC_STUDY, CELL_STUDY)))
drawn_tab[, nuclei := study %in% NUC_STUDY]
IA <- fromJSON(A[["islet_roster"]])$adult_islet
EXCL <- names(IA$excluded_units)
drop <- drawn_tab[, grepl("^Adult islet", grp) & study %in% IA$datasets & unit %in% EXCL]
stopifnot(sum(drop) == 3L, setequal(drawn_tab$grp[drop], c("Adult islet beta", "Adult islet alpha", "Adult islet delta")))   # one donor, three cell types
drawn_tab <- drawn_tab[!drop]
con <- file(A[["drawn"]], "wb"); write.table(drawn_tab, con, sep = "\t", quote = TRUE, row.names = FALSE, eol = "\n"); close(con)
cat("STAC tables written\n")
