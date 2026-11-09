# Figure S8B (cell-type column): cell-type signature enrichment of the HIO three-way gene patterns of Figure S8A.
# Lists: the six pattern x direction gene lists (stem_up/down, gut_up/down, hio_up/down). Test: one-sided hypergeometric
# (clusterProfiler::enricher) against MSigDB C8 cell-type signature gene sets (MSigDB 2026.1 via msigdbr, read from the
# release zip); Benjamini-Hochberg within list; gene sets of 10-500 genes; called at adjusted p < 0.05. Background: the
# expressed genes of Figure S8A.
# Inputs:  genes_all    genes_all.tsv.gz from s4_hio_three_way_patterns.R
#          msigdb_zip   msigdb.2026.1.zip (the msigdbr 26.1 release archive)
# Outputs: all          celltype_all_terms.tsv (every tested term)
#          called       celltype_called.tsv (terms at adjusted p < 0.05)
#   Rscript s4_hio_celltype_signatures.R genes_all=... msigdb_zip=... all=celltype_all_terms.tsv called=celltype_called.tsv
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({library(data.table); library(clusterProfiler); library(msigdbr)})

# msigdbr 26.1 reads its gene sets from the release zip unpacked in tools::R_user_dir("msigdbr", "cache") (it looks for
# the zip there and the .rds files it holds); set that up in a temporary folder from the given zip, so nothing is downloaded.
msigdb_sets <- function(zip, collection) {
  root <- tempfile("msigdbr_cache_"); Sys.setenv(R_USER_CACHE_DIR = root)
  cache <- tools::R_user_dir("msigdbr", "cache"); dir.create(cache, recursive = TRUE)
  stopifnot(unname(tools::md5sum(zip)) == "512ba99c6827141a9d471972b812d4ac")   # the release MD5 built into msigdbr 26.1
  file.copy(zip, file.path(cache, "msigdb.2026.1.zip")); utils::unzip(zip, exdir = cache)
  as.data.frame(msigdbr(species = "Homo sapiens", collection = collection))[, c("gs_name", "gene_symbol")]
}

G <- as.data.table(read.delim(gzfile(A[["genes_all"]]), stringsAsFactors = FALSE))
U <- G[expressed == TRUE, gene]
tag <- c(STEM_SHARED = "stem", GUT_SHARED = "gut", HIO_ONLY = "hio")
G[, list := fifelse(pattern == "none", NA_character_, paste0(tag[pattern], "_", fifelse(grepl("^higher", direction), "up", "down")))]
LISTS <- split(G[!is.na(list), gene], G[!is.na(list), list])
stopifnot(length(LISTS) == 6, all(unlist(LISTS) %in% U), length(U) == 11496L)

C8 <- msigdb_sets(A[["msigdb_zip"]], "C8")
stopifnot(length(unique(C8$gs_name)) == 866L)   # C8 sets in MSigDB 2026.1
P <- list(pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
RES <- rbindlist(lapply(names(LISTS), function(n) {
  r <- do.call(enricher, c(list(gene = LISTS[[n]], universe = U, TERM2GENE = C8), P)); d <- if (is.null(r)) NULL else as.data.frame(r)
  if (is.null(d) || !nrow(d)) return(NULL)
  data.table(list = n, list_genes = length(LISTS[[n]]), collection = "C8_CELL_TYPE", term_id = d$ID, term = d$Description,
             genes_in_list_and_term = d$Count, gene_ratio = d$GeneRatio, bg_ratio = d$BgRatio, p = d$pvalue, p_adj = d$p.adjust, genes = d$geneID) }))
RES[, called := p_adj < 0.05]; setorder(RES, list, p)
fwrite(RES, A[["all"]], sep = "\t"); fwrite(RES[called == TRUE], A[["called"]], sep = "\t")
