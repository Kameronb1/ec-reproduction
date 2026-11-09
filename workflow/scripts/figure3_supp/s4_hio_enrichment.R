# Figure S8B (pathway column): pathway enrichment of the HIO three-way gene patterns of Figure S8A.
# Lists: the six pattern x direction gene lists (stem_up/down, gut_up/down, hio_up/down) and the STEM_SHARED lists
# restricted to genes marked "EC only". Test: one-sided hypergeometric (clusterProfiler) against GO biological process
# (org.Hs.eg.db), Reactome (ReactomePA) and MSigDB Hallmark (MSigDB 2026.1 via msigdbr, read from the release zip);
# Benjamini-Hochberg within list x collection; gene sets of 10-500 genes; called at adjusted p < 0.05. Background: the
# expressed genes of Figure S8A.
# Inputs:  genes_all    genes_all.tsv.gz from s4_hio_three_way_patterns.R
#          msigdb_zip   msigdb.2026.1.zip (the msigdbr 26.1 release archive)
# Outputs: all          enrichment_all_terms.tsv (every tested term)
#          called       enrichment_called.tsv (terms at adjusted p < 0.05)
#   Rscript s4_hio_enrichment.R genes_all=... msigdb_zip=... all=enrichment_all_terms.tsv called=enrichment_called.tsv
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA); library(msigdbr)})

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
LISTS$stem_up_EC_only <- G[list == "stem_up" & other_endocrine == "EC only", gene]
LISTS$stem_down_EC_only <- G[list == "stem_down" & other_endocrine == "EC only", gene]
stopifnot(all(unlist(LISTS) %in% U), length(U) == 11496L)

ENTREZ <- suppressMessages(bitr(U, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)); ENTREZ <- ENTREZ[!duplicated(ENTREZ$SYMBOL), ]
HALL <- msigdb_sets(A[["msigdb_zip"]], "H")
P <- list(pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
run <- function(genes) {
  e <- ENTREZ$ENTREZID[ENTREZ$SYMBOL %in% genes]
  r <- list(GO_BP = do.call(enrichGO, c(list(gene = genes, universe = U, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP"), P)),
            REACTOME = do.call(enrichPathway, c(list(gene = e, universe = ENTREZ$ENTREZID, organism = "human", readable = TRUE), P)),
            HALLMARK = do.call(enricher, c(list(gene = genes, universe = U, TERM2GENE = HALL), P)))
  rbindlist(lapply(names(r), function(k) { d <- if (is.null(r[[k]])) NULL else as.data.frame(r[[k]])
    if (is.null(d) || !nrow(d)) return(NULL)
    data.table(collection = k, term_id = d$ID, term = d$Description, genes_in_list_and_term = d$Count, gene_ratio = d$GeneRatio,
               bg_ratio = d$BgRatio, p = d$pvalue, p_adj = d$p.adjust, genes = d$geneID) }))
}
RES <- rbindlist(lapply(names(LISTS), function(n) { r <- run(LISTS[[n]]); if (!nrow(r)) NULL else cbind(list = n, list_genes = length(LISTS[[n]]), r) }))
RES[, called := p_adj < 0.05]; setorder(RES, list, collection, p)
fwrite(RES, A[["all"]], sep = "\t"); fwrite(RES[called == TRUE], A[["called"]], sep = "\t")
