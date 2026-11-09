# Figure S13B: cell-type signature enrichment of the Figure 3 gene groups with and without their mixed genes.
# Lists: the nine lists of Figure S13A (lists.tsv from s8_mixed_gene_enrichment.R): G1, G2, G3 alone, each with the mixed
# genes it would join, and the mixed lists alone. Test: one-sided hypergeometric (clusterProfiler::enricher) against
# MSigDB C8 cell-type signature gene sets (MSigDB 2026.1 via msigdbr, read from the release zip); Benjamini-Hochberg
# within list; gene sets of 10-500 genes; called at adjusted p < 0.05; background the 14,974 genes the donor-level
# Figure 3D fit tested. Each called term gets a class by fixed name rules (intestinal EC or enteroendocrine; neuroendocrine,
# other organs; glia or neural progenitor; neuron; other). Each group is compared with itself plus its mixed genes
# (called terms alone / with mixed / kept / new / lost).
# Inputs:  lists        lists.tsv (list, gene) from s8_mixed_gene_enrichment.R
#          universe     universe_paneld_tested.txt (background genes)
#          msigdb_zip   msigdb.2026.1.zip (the msigdbr 26.1 release archive)
# Outputs: all          celltype_all_terms.tsv (every tested term of every list, with term class)
#          called       celltype_called.tsv (terms at adjusted p < 0.05)
#          compare      group_vs_group_with_mixed.tsv (alone / with_mixed / kept / new / lost per group)
#          by_class     called_terms_by_class.tsv (called terms per list and term class)
#   Rscript s8_mixed_gene_celltypes.R lists=... universe=... msigdb_zip=... all=... called=... compare=... by_class=...
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(data.table); library(clusterProfiler); library(msigdbr) })

# msigdbr 26.1 reads its gene sets from the release zip unpacked in tools::R_user_dir("msigdbr", "cache") (it looks for
# the zip there and the .rds files it holds); set that up in a temporary folder from the given zip, so nothing is downloaded.
msigdb_sets <- function(zip, collection) {
  root <- tempfile("msigdbr_cache_"); Sys.setenv(R_USER_CACHE_DIR = root)
  cache <- tools::R_user_dir("msigdbr", "cache"); dir.create(cache, recursive = TRUE)
  stopifnot(unname(tools::md5sum(zip)) == "512ba99c6827141a9d471972b812d4ac")   # the release MD5 built into msigdbr 26.1
  file.copy(zip, file.path(cache, "msigdb.2026.1.zip")); utils::unzip(zip, exdir = cache)
  as.data.frame(msigdbr(species = "Homo sapiens", collection = collection))[, c("gs_name", "gene_symbol")]
}

LT <- fread(A[["lists"]]); LISTS <- split(LT$gene, LT$list)
want <- c(G1 = 28L, G1_PLUS_MIXED_G1_LIKE = 89L, G2 = 34L, G2_PLUS_MIXED_G2_LIKE = 69L, G3 = 84L, G3_PLUS_MIXED_G3_IF_NOT_SHARED = 100L,
          MIXED_G1_LIKE = 61L, MIXED_G2_LIKE = 35L, MIXED_G3_IF_NOT_SHARED = 16L)
if (!setequal(names(LISTS), names(want)) || any(lengths(LISTS)[names(want)] != want)) stop("list sizes")
U <- readLines(A[["universe"]]); U <- U[nzchar(U)]
if (length(U) != 14974L || !all(unlist(LISTS) %in% U)) stop("background")
C8 <- msigdb_sets(A[["msigdb_zip"]], "C8")
P <- list(pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
RES <- rbindlist(lapply(names(want), function(n) {
  r <- do.call(enricher, c(list(gene = LISTS[[n]], universe = U, TERM2GENE = C8), P)); d <- if (is.null(r)) NULL else as.data.frame(r)
  if (is.null(d) || !nrow(d)) return(NULL)
  data.table(list = n, list_genes = length(LISTS[[n]]), term = d$Description, genes_in_list_and_term = d$Count, p = d$pvalue, p_adj = d$p.adjust, genes = d$geneID) }))
RES[, called := p_adj < 0.05]

# term classes (first matching rule wins)
GUT_EC <- c("BUSSLINGER_DUODENAL_EC_CELLS", "GAO_LARGE_INTESTINE_ADULT_CA_ENTEROENDOCRINE_CELLS", "DESCARTES_FETAL_INTESTINE_CHROMAFFIN_CELLS")
RULES <- list("Gut EC or enteroendocrine" = function(t) t %in% GUT_EC,
  "Neuroendocrine, other organs" = function(t) grepl("NEUROENDOCRINE|_NE_PRECURSOR", t),
  "Glia or neural progenitor" = function(t) grepl("_OPC$|_HOPC$|_HRGL|OLIGODENDROCYTES|MICROGLIA|_HNPROG$", t),
  "Neuron" = function(t) grepl("NEURON|SYMPATHOBLASTS|MANNO_MIDBRAIN_NEUROTYPES_H(DA|DA2|SERT|OMTN|NBML1|NBML5)$|PURKINJE", t),
  "Other" = function(t) rep(TRUE, length(t)))
cls <- function(t) { out <- rep(NA_character_, length(t)); for (k in names(RULES)) { hit <- is.na(out) & RULES[[k]](t); out[hit] <- k }; out }
RES[, term_class := factor(cls(term), names(RULES))]
setorder(RES, list, p)
fwrite(RES, A[["all"]], sep = "\t"); fwrite(RES[called == TRUE], A[["called"]], sep = "\t")

# kept / new / lost per group, and class counts
PAIR <- c(G1 = "G1_PLUS_MIXED_G1_LIKE", G2 = "G2_PLUS_MIXED_G2_LIKE", G3 = "G3_PLUS_MIXED_G3_IF_NOT_SHARED")
CMP <- rbindlist(lapply(names(PAIR), function(g) { a <- RES[list == g & called == TRUE, term]; b <- RES[list == PAIR[[g]] & called == TRUE, term]
  data.table(group = g, alone = length(a), with_mixed = length(b), kept = length(intersect(a, b)), new = length(setdiff(b, a)), lost = length(setdiff(a, b))) }))
CLS <- dcast(RES[called == TRUE, .N, by = .(list, term_class)], list ~ term_class, value.var = "N", fill = 0)
fwrite(CMP, A[["compare"]], sep = "\t"); fwrite(CLS, A[["by_class"]], sep = "\t")
