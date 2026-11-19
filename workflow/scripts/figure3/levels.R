# Figure 3A and 3B: expression of the three gene groups in every population, from per-sample pseudobulk sums.
#
# Level of a gene in a profile = log2(1 + CPM) of the profile's summed raw counts, CPM over the profile's full-cell library.
#   3A  one profile per sample and population (samples with at least 10 cells of the population); per profile and group,
#       the median level over the group's genes. Also the median over the group's genes of (level - primary intestinal EC
#       reference), the reference being the per-gene mean of the three intestinal studies' pooled levels (Elmentaite 2020
#       and 2021 averaged first).
#   3B  the 28 G1 genes, counts pooled over all samples of a study's population (intestinal EC, HIO graft EC, SC-SICs).
# The seven Elmentaite 2020 donors with Crohn's disease and the adult islet donor listed as excluded in the adult islet
# roster are left out.
#
# Inputs:  groups (gene groups), {gut,hio,sc,islet}_sums and _groups (pseudobulk sums and profile tables),
#          gut_roster, islet_roster
# Outputs: medians (3A), g1_levels (3B), g1_counts (samples and cells behind each 3B column)
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

GG <- as.data.table(read.delim(gzfile(A[["groups"]]), check.names = FALSE, stringsAsFactors = FALSE))
GRPS <- c(G1 = "G1_EC_GENES_SC_SIC_LACK", G2 = "G2_EC_GENES_SC_SIC_HAVE", G3 = "G3_SC_SIC_OWN_PROGRAM")
GG[, g := names(GRPS)[match(group, GRPS)]]
stopifnot(!anyNA(GG$g), identical(as.vector(table(GG$g)[names(GRPS)]), c(28L, 34L, 84L)), !anyDuplicated(GG$gene))
GENES <- GG$gene

# pseudobulk sums: only the 146 group genes' rows are kept
read_src <- function(k) {
  g <- fread(A[[paste0(k, "_groups")]], colClasses = list(character = "unit"))
  s <- as.data.table(read.delim(gzfile(A[[paste0(k, "_sums")]]), check.names = FALSE, stringsAsFactors = FALSE))
  stopifnot(identical(names(s)[-1], g$group), !anyDuplicated(s$gene))
  miss <- setdiff(GENES, s$gene); if (length(miss)) stop(k, ": group genes missing from the sums: ", paste(miss, collapse = ", "))
  m <- as.matrix(s[match(GENES, gene), -1]); rownames(m) <- GENES
  list(g = g, m = m)
}
SRCS <- c(gut = "gut", hio = "hio", sc = "sc", islet = "islet")
D <- lapply(SRCS, read_src)

crohns <- fromJSON(A[["gut_roster"]])$excluded_units
IR <- fromJSON(A[["islet_roster"]])$adult_islet
stopifnot(length(crohns) == 7L, identical(names(IR$excluded_units), "HPAP034"))

# populations
POP <- rbindlist(lapply(names(SRCS), function(k) cbind(D[[k]]$g, src = k)), fill = TRUE)
POP <- POP[!(src == "gut" & family_id == "SF_ELMENTAITE2020_DEVGUT" & unit %in% crohns)]
POP <- POP[!(src == "islet" & family_id %in% IR$datasets & unit %in% names(IR$excluded_units))]
POP[, pop := fcase(src == "gut" & class == "POSITIVE", "Primary intestinal EC", src == "gut" & class == "NEGATIVE", "Other primary EEC",
                   src == "hio" & class == "EC", "HIO graft EC",
                   src == "sc" & class == "SC_SIC", "SC-SIC", src == "sc" & class == "SC_BETA", "SC-beta",
                   src == "sc" & class == "SC_ALPHA", "SC-alpha", src == "sc" & class == "SC_DELTA", "SC-delta",
                   src == "islet" & class %in% c("Beta", "beta_major", "beta_minor"), "Adult islet beta",
                   src == "islet" & class %in% c("Alpha", "alpha", "Proliferating_Alpha"), "Adult islet alpha", default = NA_character_)]
POP <- POP[!is.na(pop)]
umi_of <- function(k, groups) rowSums(D[[k]]$m[, groups, drop = FALSE])
# one profile per sample and population (beta_major + beta_minor and Alpha + Proliferating_Alpha summed per donor)
S <- POP[, .(n_cells = sum(n_cells), lib = sum(lib), groups = list(group)), by = .(pop, src, study = family_id, unit)]
S <- S[n_cells >= 10L]
S[, sample_id := paste(study, unit, pop, sep = "|")]
LV <- sapply(seq_len(nrow(S)), function(i) log2(1 + umi_of(S$src[i], S$groups[[i]]) / S$lib[i] * 1e6))
colnames(LV) <- S$sample_id
S[, nuclei := study %in% c("GSE199636", "SF_HICKEY2023_ADULT_EEC")]          # single-nucleus studies
DEFINING <- c("Primary intestinal EC", "Other primary EEC", "SC-SIC", "SC-beta", "SC-alpha", "SC-delta")
S[, used_to_define := pop %in% DEFINING]

# primary intestinal EC reference per gene: mean of (Elmentaite 2020 + 2021 mean), Hickey 2023, Egozi 2021
pooled <- function(k, groups, lib) log2(1 + umi_of(k, groups) / sum(lib) * 1e6)
fam_pool <- function(fam) { x <- POP[pop == "Primary intestinal EC" & family_id == fam]; pooled("gut", x$group, x$lib) }
REF <- ((fam_pool("SF_ELMENTAITE2021_GUT_ATLAS") + fam_pool("SF_ELMENTAITE2020_DEVGUT")) / 2 +
          fam_pool("SF_HICKEY2023_ADULT_EEC") + fam_pool("SF_EGOZI2021_DEVELOPMENTAL_EEC")) / 3

# 3A: per profile and group, the median over the group's genes
long <- data.table(sample_id = rep(colnames(LV), each = nrow(LV)), gene = rep(rownames(LV), ncol(LV)), level = as.vector(LV))
long <- merge(S[, .(sample_id, pop, study, unit, n_cells, lib, nuclei, used_to_define)], long, by = "sample_id")
long[, g := GG$g[match(gene, GG$gene)]]
long[, gap := level - REF[gene]]
med <- long[, .(genes = .N, median_level = median(level), median_gap_to_primary_ec = median(gap)),
            by = .(sample_id, pop, study, unit, n_cells, nuclei, used_to_define, g)]
setorder(med, pop, study, unit, g)
fwrite(med, A[["medians"]], sep = "\t")

# 3B: the 28 G1 genes, counts summed over a study's samples divided by the summed full-cell library
g1 <- GG[g == "G1"][order(-abs(paneld_t))]$gene
col_pool <- function(src_k, popname, fams) { x <- POP[src == src_k & pop == popname & family_id %in% fams]
  stopifnot(nrow(x) > 0); pooled(src_k, x$group, x$lib)[g1] }
SC_STUDY <- c(GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022", GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024",
              GSE185036 = "Du 2022", GSE199636 = "Augsornworawat 2023")
B <- data.table(gene = g1,
  `Elmentaite 2020 + 2021` = col_pool("gut", "Primary intestinal EC", c("SF_ELMENTAITE2020_DEVGUT", "SF_ELMENTAITE2021_GUT_ATLAS")),
  `Hickey 2023` = col_pool("gut", "Primary intestinal EC", "SF_HICKEY2023_ADULT_EEC"),
  `Egozi 2021` = col_pool("gut", "Primary intestinal EC", "SF_EGOZI2021_DEVELOPMENTAL_EEC"),
  `McCauley 2023` = col_pool("hio", "HIO graft EC", "SF_GSE214852_HIO_GRAFT"))
for (k in names(SC_STUDY)) B[, (SC_STUDY[[k]]) := col_pool("sc", "SC-SIC", k)]
B[, row := .I]
fwrite(B, A[["g1_levels"]], sep = "\t")

# samples and cells behind each 3B column
BN <- rbindlist(list(
  POP[pop == "Primary intestinal EC", .(column = fcase(family_id %in% c("SF_ELMENTAITE2020_DEVGUT", "SF_ELMENTAITE2021_GUT_ATLAS"), "Elmentaite 2020 + 2021",
        family_id == "SF_HICKEY2023_ADULT_EEC", "Hickey 2023", default = "Egozi 2021"), unit, n_cells)],
  POP[pop == "HIO graft EC", .(column = "McCauley 2023", unit, n_cells)],
  POP[pop == "SC-SIC", .(column = SC_STUDY[family_id], unit, n_cells)]))
fwrite(BN[, .(samples = uniqueN(unit), cells = sum(n_cells)), by = column], A[["g1_counts"]], sep = "\t")
