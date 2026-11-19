# Figure S4A and S4B: expression of the 34 G2 and 84 G3 genes in each study, log2(1 + CPM) of raw counts summed over the
# study's samples of the population, CPM over the summed full-cell library. Columns: primary intestinal EC of three
# studies (Elmentaite 2020 and 2021 pooled), McCauley 2023 HIO graft EC, SC-SICs of six SC-islet studies. The seven
# Elmentaite 2020 donors with Crohn's disease are left out.
# Gene order as drawn: G2 by the SC-SIC minus intestinal EC log fold change of the donor-level fit (largest first), split
# into the genes called higher in SC-SICs and those without a call; G3 by the size of the study-level t statistic.
#
# Inputs:  groups (gene groups), {gut,hio,sc}_sums and _groups (pseudobulk sums and profile tables), gut_roster
# Outputs: levels (gene x study), counts (samples and cells per column), order (gene order and blocks)
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

GG <- as.data.table(read.delim(gzfile(A[["groups"]]), check.names = FALSE, stringsAsFactors = FALSE))
GRPS <- c(G1 = "G1_EC_GENES_SC_SIC_LACK", G2 = "G2_EC_GENES_SC_SIC_HAVE", G3 = "G3_SC_SIC_OWN_PROGRAM")
GG[, g := names(GRPS)[match(group, GRPS)]]
GENES <- GG[g %in% c("G2", "G3")]$gene
stopifnot(length(GENES) == 118L)

read_src <- function(k) {
  g <- fread(A[[paste0(k, "_groups")]], colClasses = list(character = "unit"))
  s <- as.data.table(read.delim(gzfile(A[[paste0(k, "_sums")]]), check.names = FALSE, stringsAsFactors = FALSE))
  stopifnot(identical(names(s)[-1], g$group), !anyDuplicated(s$gene))
  miss <- setdiff(GENES, s$gene); if (length(miss)) stop(k, ": genes missing from the sums: ", paste(miss, collapse = ", "))
  m <- as.matrix(s[match(GENES, gene), -1]); rownames(m) <- GENES
  list(g = g, m = m)
}
SRCS <- c(gut = "gut", hio = "hio", sc = "sc")
D <- lapply(SRCS, read_src)
crohns <- fromJSON(A[["gut_roster"]])$excluded_units
stopifnot(length(crohns) == 7L)
POP <- rbindlist(lapply(names(SRCS), function(k) cbind(D[[k]]$g, src = k)), fill = TRUE)
POP <- POP[!(src == "gut" & family_id == "SF_ELMENTAITE2020_DEVGUT" & unit %in% crohns)]
POP[, pop := fcase(src == "gut" & class == "POSITIVE", "Primary intestinal EC", src == "hio" & class == "EC", "HIO graft EC",
                   src == "sc" & class == "SC_SIC", "SC-SIC", default = NA_character_)]
POP <- POP[!is.na(pop)]
pooled <- function(k, groups, lib) log2(1 + rowSums(D[[k]]$m[, groups, drop = FALSE]) / sum(lib) * 1e6)
col_pool <- function(src_k, popname, fams) { x <- POP[src == src_k & pop == popname & family_id %in% fams]; stopifnot(nrow(x) > 0); pooled(src_k, x$group, x$lib) }
L <- data.table(gene = GENES, g = GG$g[match(GENES, GG$gene)],
  `Elmentaite 2020 + 2021` = col_pool("gut", "Primary intestinal EC", c("SF_ELMENTAITE2020_DEVGUT", "SF_ELMENTAITE2021_GUT_ATLAS")),
  `Egozi 2021` = col_pool("gut", "Primary intestinal EC", "SF_EGOZI2021_DEVELOPMENTAL_EEC"),
  `Hickey 2023` = col_pool("gut", "Primary intestinal EC", "SF_HICKEY2023_ADULT_EEC"),
  `McCauley 2023` = col_pool("hio", "HIO graft EC", "SF_GSE214852_HIO_GRAFT"),
  `Augsornworawat 2020` = col_pool("sc", "SC-SIC", "GSE151117"), `Balboa 2022` = col_pool("sc", "SC-SIC", "GSE167880"),
  `Du 2022` = col_pool("sc", "SC-SIC", "GSE185036"), `Augsornworawat 2023` = col_pool("sc", "SC-SIC", "GSE199636"),
  `Zhu 2023` = col_pool("sc", "SC-SIC", "GSE202497"), `Siehler 2024` = col_pool("sc", "SC-SIC", "GSE226346"))
fwrite(L, A[["levels"]], sep = "\t")

FAM <- c(SF_ELMENTAITE2020_DEVGUT = "Elmentaite 2020 + 2021", SF_ELMENTAITE2021_GUT_ATLAS = "Elmentaite 2020 + 2021", SF_EGOZI2021_DEVELOPMENTAL_EEC = "Egozi 2021",
         SF_HICKEY2023_ADULT_EEC = "Hickey 2023", SF_GSE214852_HIO_GRAFT = "McCauley 2023", GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022",
         GSE185036 = "Du 2022", GSE199636 = "Augsornworawat 2023", GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024")
CN <- POP[family_id %in% names(FAM), .(samples = uniqueN(unit), cells = sum(n_cells), crohn_units = sum(unit %in% crohns & family_id == "SF_ELMENTAITE2020_DEVGUT")),
          by = .(column = FAM[family_id])]
CN <- CN[, .(samples = sum(samples), cells = sum(cells), crohn_units = sum(crohn_units)), by = column]
fwrite(CN, A[["counts"]], sep = "\t")

# gene order of the two panels
g2 <- GG[g == "G2"][order(-panel5_logFC)]
g3 <- GG[g == "G3"][order(-abs(paneld_t))]
O <- rbind(data.table(group = "G2", part = "higher in SC-SICs", gene = g2[paneld_call == "SC_SIC_HIGHER"]$gene),
           data.table(group = "G2", part = "no significant difference", gene = g2[paneld_call == "NO_CALL"]$gene),
           data.table(group = "G3", part = "all", gene = g3$gene))
stopifnot(identical(O[, .N, by = part]$N, c(9L, 25L, 84L)))
O[, position := seq_len(.N), by = group]
fwrite(O, A[["order"]], sep = "\t")
