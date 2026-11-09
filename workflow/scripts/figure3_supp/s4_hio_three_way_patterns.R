# Figure S8A: three-way comparison of EC gene expression between HIO graft EC (H), SC-SIC (S) and primary intestinal EC (P).
# For every gene shared by the three pseudobulk inputs, a unit's value is log2 CPM (edgeR TMM over all units together,
# prior count 1), restricted to protein-coding genes on chr1-22 and X (org.Hs.eg.db) and corrected for a per-unit trend
# with genomic gene span (loess of the unit's log2 CPM minus the all-unit gene mean on log10 span, fitted on genes with
# mean log2 CPM >= 1, then subtracted). A study's value is the mean of its units; a group's value is the mean of its
# studies. Groups: P = primary intestinal EC (Egozi 2021, Elmentaite 2020+2021 as one study, Hickey 2023; the seven Elmentaite
# 2020 child Crohn's donors removed), S = SC-SIC (six studies, ten samples), H = HIO graft EC (McCauley 2023
# "Endocrine_1" primary study label, two organoids). Units with fewer than 10 cells are dropped.
# Each expressed gene (max group mean log2 CPM >= 2) gets one pattern, set by its closest pair of groups (T = 1 log2):
#   STEM_SHARED  H and S both differ from P (and from fetal intestinal EC) in the same direction by >= T; both organoids
#                and >= 5 of 6 SC-SIC studies beyond every intestinal study mean
#   GUT_SHARED   S differs from both P and H by >= T in the same direction; both organoids closer to P than to S;
#                >= 5 of 6 SC-SIC studies beyond every intestinal study mean
#   HIO_ONLY     H differs from both P and S by >= T in the same direction; both organoids beyond every intestinal and
#                SC-SIC study mean
# Also: the same differences in the non-EC endocrine cells ("also in other endocrine" / "EC only"), and a sensitivity
# run without the single-nucleus SC-SIC study GSE199636 (in_sensitivity). No p-values (H is two organoids of one line);
# genes are ranked by effect size (score).
# Inputs:  gut_sums, gut_groups   primary intestinal endocrine pseudobulk raw-count sums per unit x class, and their group table
#          sc_sums, sc_groups     SC-islet endocrine pseudobulk sums and group table
#          hio_sums, hio_groups   HIO (McCauley 2023) endocrine pseudobulk sums and group table
#          roster                 intestinal reference donor list (JSON; excluded_units = the Crohn's donors)
# Output:  out                    genes_all.tsv.gz, one row per tested gene with group means, differences and pattern
#   Rscript s4_hio_three_way_patterns.R gut_sums=... gut_groups=... sc_sums=... sc_groups=... hio_sums=... hio_groups=... roster=... out=genes_all.tsv.gz
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({library(data.table); library(edgeR); library(org.Hs.eg.db); library(jsonlite)})
T_LOG2 <- 1; FLOOR <- 10L; MIN_S_STUDIES <- 5L
E20 <- "SF_ELMENTAITE2020_DEVGUT"; E21 <- "SF_ELMENTAITE2021_GUT_ATLAS"; EGOZI <- "SF_EGOZI2021_DEVELOPMENTAL_EEC"; HICKEY <- "SF_HICKEY2023_ADULT_EEC"
SNUC <- "GSE199636"

rd <- function(f) { m <- as.data.table(read.delim(gzfile(f), check.names = FALSE)); g <- m$gene; m <- as.matrix(m[, -1]); rownames(m) <- g; stopifnot(!anyDuplicated(g)); m }
gs <- rd(A[["gut_sums"]]); gg <- fread(A[["gut_groups"]])
ss <- rd(A[["sc_sums"]]); sg <- fread(A[["sc_groups"]])
hs <- rd(A[["hio_sums"]]); hg <- fread(A[["hio_groups"]])
stopifnot(identical(colnames(gs), gg$group), identical(colnames(ss), sg$group), identical(colnames(hs), hg$group))
for (m in list(gs, ss, hs)) stopifnot(all(m >= 0), all(m == round(m)))
genes <- Reduce(intersect, list(rownames(gs), rownames(ss), rownames(hs)))
stopifnot(length(genes) == nrow(gs), nrow(gs) == 21520)   # the intestinal input holds only the 21,520 shared genes

# Units: system, study, unit, class (EC / OTHER), age (intestinal units only).
crohns <- fromJSON(A[["roster"]])$excluded_units; stopifnot(length(crohns) == 7)
gg <- gg[!(family_id == E20 & unit %in% crohns)]
gg[, `:=`(system = "P", study = fifelse(family_id %in% c(E20, E21), "Elmentaite 2020+2021", fifelse(family_id == EGOZI, "Egozi 2021", "Hickey 2023")),
          cls = fifelse(class == "POSITIVE", "EC", "OTHER"),
          age = fcase(family_id == HICKEY, "adult",
                      family_id == EGOZI, fifelse(grepl("_(21|23)w$", unit), "fetal", "neonatal"),
                      family_id == E21, fifelse(startsWith(unit, "F"), "fetal", "adult"),
                      family_id == E20, fifelse(grepl("^[0-9]", unit), "fetal", "child")))]
stopifnot(!anyNA(gg$age))
# SC-islet other endocrine: alpha + beta + delta summed per sample.
so <- sg[class %in% c("SC_ALPHA", "SC_BETA", "SC_DELTA")]
so_sum <- sapply(split(so$group, so$unit), function(g) rowSums(ss[, g, drop = FALSE]))
so_meta <- so[, .(n_cells = sum(n_cells), family_id = family_id[1]), by = unit]; so_sum <- so_sum[, so_meta$unit, drop = FALSE]
colnames(so_sum) <- paste0(so_meta$family_id, "|", so_meta$unit, "|SC_OTHER_ENDOCRINE")
sic <- sg[class == "SC_SIC"]
stopifnot(nrow(sic) == 10, uniqueN(sic$family_id) == 6)
counts <- cbind(gs[genes, gg$group], ss[genes, sic$group], so_sum[genes, ], hs[genes, hg$group])
meta <- rbind(gg[, .(group, system, study, unit, cls, age, n_cells)],
              sic[, .(group, system = "S", study = family_id, unit, cls = "EC", age = NA_character_, n_cells)],
              data.table(group = colnames(so_sum), system = "S", study = so_meta$family_id, unit = so_meta$unit, cls = "OTHER", age = NA_character_, n_cells = so_meta$n_cells),
              hg[, .(group, system = "H", study = "McCauley 2023", unit, cls = fifelse(class == "EC", "EC", "OTHER"), age = NA_character_, n_cells)])
stopifnot(identical(meta$group, colnames(counts)))
stopifnot(identical(hg[class == "EC"]$n_cells, c(180L, 142L)), identical(hg[class == "OTHER_EEC"]$n_cells, c(200L, 76L)))
keep <- meta$n_cells >= FLOOR; meta <- meta[keep]; counts <- counts[, keep]

# log2 CPM per unit (TMM over all units together), then study means, then group means.
y <- normLibSizes(DGEList(counts)); lc <- cpm(y, log = TRUE, prior.count = 1)
# Gene annotation: type and genomic span (largest end minus smallest start on the gene's main chromosome).
sym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = rownames(lc), keytype = "SYMBOL", columns = c("ENTREZID", "GENETYPE")))
sym <- as.data.table(sym)[!is.na(ENTREZID)][, .SD[1], by = SYMBOL]
loc <- merge(as.data.table(toTable(org.Hs.egCHRLOC)), as.data.table(toTable(org.Hs.egCHRLOCEND)), by = c("gene_id", "Chromosome"), allow.cartesian = TRUE)
loc <- loc[Chromosome %in% c(1:22, "X", "Y"), .(span = max(abs(end_location)) - min(abs(start_location)), chrom = Chromosome[1]), by = gene_id]
ann <- merge(sym, loc, by.x = "ENTREZID", by.y = "gene_id")
ann <- ann[GENETYPE == "protein-coding" & chrom %in% c(1:22, "X") & span > 0]
kept_genes <- intersect(rownames(lc), ann$SYMBOL); span <- log10(ann$span[match(kept_genes, ann$SYMBOL)])
lc <- lc[kept_genes, ]
# Per-unit span trend removed (each unit against the all-unit gene mean).
gm <- rowMeans(lc); fitg <- gm >= 1
trend <- sapply(seq_len(ncol(lc)), function(j) {
  r <- lc[, j] - gm; f <- loess(r ~ s, data = data.frame(r = r[fitg], s = span[fitg]), span = .5, degree = 1)
  predict(f, data.frame(s = pmin(pmax(span, min(span[fitg])), max(span[fitg]))))
})
lc <- lc - trend
study_means <- function(cols) { m <- meta[cols]; sapply(split(which(cols), m$study), function(i) rowMeans(lc[, i, drop = FALSE])) }
grp <- function(sys, cl, extra = TRUE) { cols <- meta$system == sys & meta$cls == cl & extra; study_means(cols) }
run <- function(drop_snuc) {
  sx <- if (drop_snuc) meta$study != SNUC else TRUE
  P <- grp("P", "EC"); Pf <- grp("P", "EC", meta$age %in% "fetal"); S <- grp("S", "EC", sx)
  Hu <- lc[, meta$system == "H" & meta$cls == "EC"]; colnames(Hu) <- meta[system == "H" & cls == "EC"]$unit
  Po <- grp("P", "OTHER"); So <- grp("S", "OTHER", sx); Ho <- lc[, meta$system == "H" & meta$cls == "OTHER"]
  d <- data.table(gene = rownames(lc), P = rowMeans(P), P_fetal = rowMeans(Pf), S = rowMeans(S), H = rowMeans(Hu),
                  H_wt1 = Hu[, "wt1"], H_wt2 = Hu[, "wt2"], P_min = apply(P, 1, min), P_max = apply(P, 1, max),
                  S_min = apply(S, 1, min), S_max = apply(S, 1, max), n_S_studies = ncol(S),
                  S_up_n = rowSums(S > apply(P, 1, max)), S_down_n = rowSums(S < apply(P, 1, min)),
                  P_other = rowMeans(Po), S_other = rowMeans(So), H_other = rowMeans(Ho))
  d[, `:=`(HP = H - P, SP = S - P, HS = H - S, HPf = H - P_fetal, SPf = S - P_fetal,
           HP_other = H_other - P_other, SP_other = S_other - P_other, HS_other = H_other - S_other)]
  # expressed: group mean log2 CPM >= 2 (about 3 CPM) in at least one EC group
  d[, expressed := pmax(P, S, H) >= 2]
  need_S <- min(MIN_S_STUDIES, ncol(S) - 1L)
  # the odd group out: the pattern is set by which two groups are closest, so the three patterns cannot overlap
  d[, closest := c("HS", "HP", "SP")[max.col(-cbind(abs(HS), abs(HP), abs(SP)), ties.method = "first")]]
  d[, stem := expressed & closest == "HS" & sign(HP) == sign(SP) & abs(HP) >= T_LOG2 & abs(SP) >= T_LOG2 &
              sign(HPf) == sign(HP) & sign(SPf) == sign(HP) & abs(HPf) >= T_LOG2 & abs(SPf) >= T_LOG2 &
              ifelse(HP > 0, pmin(H_wt1, H_wt2) > P_max & S_up_n >= need_S, pmax(H_wt1, H_wt2) < P_min & S_down_n >= need_S)]
  d[, gut_shared := expressed & closest == "HP" & sign(SP) == -sign(HS) & abs(SP) >= T_LOG2 & abs(HS) >= T_LOG2 &
                    abs(H_wt1 - P) < abs(H_wt1 - S) & abs(H_wt2 - P) < abs(H_wt2 - S) &
                    ifelse(SP > 0, S_up_n >= need_S, S_down_n >= need_S)]
  d[, hio_only := expressed & closest == "SP" & sign(HP) == sign(HS) & abs(HP) >= T_LOG2 & abs(HS) >= T_LOG2 &
                  ifelse(HP > 0, pmin(H_wt1, H_wt2) > pmax(P_max, S_max), pmax(H_wt1, H_wt2) < pmin(P_min, S_min))]
  stopifnot(!any(d$stem & d$gut_shared), !any(d$stem & d$hio_only), !any(d$gut_shared & d$hio_only))
  d[, pattern := fcase(stem, "STEM_SHARED", gut_shared, "GUT_SHARED", hio_only, "HIO_ONLY", default = "none")]
  d[, direction := fcase(pattern == "STEM_SHARED", fifelse(HP > 0, "higher in HIO and SC-SIC", "lower in HIO and SC-SIC"),
                         pattern == "GUT_SHARED", fifelse(SP > 0, "higher in SC-SIC", "lower in SC-SIC"),
                         pattern == "HIO_ONLY", fifelse(HP > 0, "higher in HIO", "lower in HIO"), default = "")]
  # is the same difference present in the non-EC endocrine cells?
  d[, other_endocrine := fcase(
    pattern == "STEM_SHARED", fifelse(sign(HP_other) == sign(HP) & sign(SP_other) == sign(HP) & abs(HP_other) >= T_LOG2 & abs(SP_other) >= T_LOG2, "also in other endocrine", "EC only"),
    pattern == "GUT_SHARED", fifelse(sign(SP_other) == sign(SP) & sign(-HS_other) == sign(SP) & abs(SP_other) >= T_LOG2 & abs(HS_other) >= T_LOG2, "also in other endocrine", "EC only"),
    pattern == "HIO_ONLY", fifelse(sign(HP_other) == sign(HP) & sign(HS_other) == sign(HP) & abs(HP_other) >= T_LOG2 & abs(HS_other) >= T_LOG2, "also in other endocrine", "EC only"),
    default = "")]
  d[, score := fcase(pattern == "STEM_SHARED", pmin(abs(HP), abs(SP)), pattern == "GUT_SHARED", pmin(abs(SP), abs(HS)),
                     pattern == "HIO_ONLY", pmin(abs(HP), abs(HS)), default = 0)]
  d[order(pattern, direction, -score)]
}
main <- run(FALSE); sens <- run(TRUE)
main[, in_sensitivity := pattern != "none" & pattern == sens$pattern[match(gene, sens$gene)]]
num <- names(main)[vapply(main, is.numeric, TRUE)]; main[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(main, A[["out"]], sep = "\t")
