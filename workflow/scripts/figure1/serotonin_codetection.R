# Figure 1D: transcripts for making and storing serotonin in the same cell.
# Per cell, is each of TPH1 (makes 5-HTP), DDC (makes serotonin) and SLC18A1 / VMAT1 (EC storage transporter) detected,
# and are all three detected together? Detection is >= 0.5 SoupX-corrected count (unrounded). The share of cells is
# reported per study and cell group:
#   Intestine    intestinal EC (primary intestine EC label; table group "Gut EC"),
#                HIO graft EC (McCauley 2023 Endocrine_1; table group "HIO EC", which the figure labels "HIO graft EC"),
#                intestinal other EEC (primary intestine "Other EEC" label; table group "Gut, other EEC")
#   SC-islet     SC-SIC, SC-beta, SC-alpha, SC-delta
#   Adult islet  adult beta, alpha, delta
# The 7 Crohn's disease donors are left out. A cell without corrected counts is left out; in adult islet these are the
# cells of the one capture SoupX could not correct (Fasolino 2022 GSM4453624). Shares before correction ("_raw") use
# the correction's own raw counts. Groups with < 30 cells are marked drawn = FALSE.
# Inputs: cells      per-study display cells (study, context, sample, published label, class, uncorrected TPH1/DDC/SLC18A1);
#                    the context "Native intestine" is the primary intestine
#         sc         SoupX raw and corrected counts, SC-islet cells
#         gut_ec     SoupX raw and corrected counts, intestinal EC cells (column corrected = SoupX applied)
#         gut_other  SoupX raw and corrected counts, primary intestine other EEC cells
#         adult      SoupX raw and corrected counts, adult islet endocrine cells
# Output: out        one row per study and group
#   Rscript serotonin_codetection.R cells=... sc=... gut_ec=... gut_other=... adult=... out=<serotonin_codetection.tsv>
args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
stopifnot(all(c("cells", "sc", "gut_ec", "gut_other", "adult", "out") %in% names(arg)))

KIT <- c("TPH1", "DDC", "SLC18A1")
CROHNS <- c("T017", "T019", "T176", "T189", "T197", "T202", "T203")
MIN_CELLS <- 30
MIN_COUNT <- 0.5                                                   # detected: corrected count >= 0.5 (raw counts: >= 1)
UNCORRECTED_CAPTURE <- "GSE148073|GSM4453624|"                     # adult islet capture SoupX could not correct

# ---- cells and groups ------------------------------------------------------------------------------------------------
cl <- read.delim(gzfile(arg[["cells"]]), check.names = FALSE, stringsAsFactors = FALSE)
cl$study[cl$study == "Elmentaite 2020 + Elmentaite 2021"] <- "Elmentaite 2020 + 2021"
stopifnot(all(cl$study[cl$sample %in% CROHNS] == "Elmentaite 2020 + 2021"))
cl <- cl[!cl$sample %in% CROHNS, ]
nat <- cl$context == "Native intestine"
cl$grp <- NA_character_
cl$grp[nat & cl$class == "EC"] <- "Gut EC"
cl$grp[cl$study == "McCauley 2023" & cl$author_label == "Endocrine_1"] <- "HIO EC"
cl$grp[nat & cl$class == "Other EEC"] <- "Gut, other EEC"
for (k in c("SC-SIC", "SC-beta", "SC-alpha", "SC-delta")) cl$grp[cl$context == "SC-islet" & cl$class == k] <- k
for (k in c("Beta", "Alpha", "Delta")) cl$grp[cl$context == "Adult islet" & cl$class == k] <- paste("Adult", tolower(k))
cl <- cl[!is.na(cl$grp), ]

# ---- corrected counts ------------------------------------------------------------------------------------------------
rd <- function(k) read.delim(gzfile(arg[[k]]), stringsAsFactors = FALSE, check.names = FALSE)
sx <- rd("sc"); ge <- rd("gut_ec"); go <- rd("gut_other"); ad <- rd("adult")
ge <- ge[ge$corrected, ]; go <- go[go$corrected, ]
cols <- paste0(rep(KIT, each = 2), c("_raw", "_soupx"))
CORR <- rbind(data.frame(cell = sx$cell, sx[, cols], check.names = FALSE),
              data.frame(cell = ge$cell, ge[, cols], check.names = FALSE),
              data.frame(cell = go$cell, go[, cols], check.names = FALSE),
              data.frame(cell = ad$cell, ad[, cols], check.names = FALSE))
stopifnot(!anyDuplicated(CORR$cell))
m <- match(cl$cell_id, CORR$cell)
cl$corrected <- !is.na(m)
stopifnot(all(cl$corrected[cl$context == "SC-islet"]), all(cl$corrected[cl$grp %in% c("Gut EC", "HIO EC")]),
          all(startsWith(cl$cell_id[cl$context == "Adult islet" & !cl$corrected], UNCORRECTED_CAPTURE)))
for (g in KIT) { cl[[paste0(g, "_raw_on")]] <- CORR[[paste0(g, "_raw")]][m] >= MIN_COUNT
  cl[[paste0(g, "_on")]] <- CORR[[paste0(g, "_soupx")]][m] >= MIN_COUNT }
cl <- cl[cl$corrected, ]
cl$kit <- cl$TPH1_on & cl$DDC_on & cl$SLC18A1_on
cl$kit_raw <- cl$TPH1_raw_on & cl$DDC_raw_on & cl$SLC18A1_raw_on

# ---- per study and group ---------------------------------------------------------------------------------------------
S <- do.call(rbind, lapply(split(cl, list(cl$study, cl$grp), drop = TRUE), function(e) data.frame(
  study = e$study[1], grp = e$grp[1], corrected = all(e$corrected), cells = nrow(e), n_all_three = sum(e$kit), all_three = mean(e$kit),
  TPH1 = mean(e$TPH1_on), DDC = mean(e$DDC_on), SLC18A1 = mean(e$SLC18A1_on),
  all_three_raw = mean(e$kit_raw), TPH1_raw = mean(e$TPH1_raw_on), DDC_raw = mean(e$DDC_raw_on), SLC18A1_raw = mean(e$SLC18A1_raw_on))))
S$drawn <- S$cells >= MIN_CELLS
stopifnot(all(S$TPH1 <= S$TPH1_raw + 1e-12 & S$DDC <= S$DDC_raw + 1e-12 & S$SLC18A1 <= S$SLC18A1_raw + 1e-12))
write.table(S, arg[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
