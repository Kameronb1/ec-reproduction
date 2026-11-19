# Figure 2E (pathways): over-representation of GO biological process, Reactome and MSigDB Hallmark terms in
# the Figure 2D gene groups.
#
# One-sided hypergeometric test (clusterProfiler) against the genes the Figure 2C fit tested; sets of 10-500
# background genes; Benjamini-Hochberg within each list x collection; called at adjusted p < 0.05. Lists:
# G1, G2, G3, G3's strict core (up in no intestinal study), and all Figure 2C calls on each side for context.
# GO from org.Hs.eg.db, Reactome from ReactomePA, Hallmark from MSigDB 2026.1 read by msigdbr.
#
# Inputs:  groups= universe= fit=  (gene_groups.py, donor_fit.R)
#          msigdb_zip=              msigdb.2026.1.zip (Zenodo record 18968178, as msigdbr 26.1.1 downloads it)
# Outputs: all_terms= called=
#
#   Rscript enrich_pathways.R groups=... universe=... fit=... msigdb_zip=... all_terms=... called=...

args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA); library(msigdbr) })

# msigdbr reads an unpacked copy of the archive from its cache folder; use a private one
cache <- file.path(tempdir(), "msigdbr_cache"); mdir <- file.path(cache, "R", "msigdbr")
dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
if (unname(tools::md5sum(arg[["msigdb_zip"]])) != "512ba99c6827141a9d471972b812d4ac") stop("msigdb.2026.1.zip has the wrong checksum")
file.copy(arg[["msigdb_zip"]], file.path(mdir, "msigdb.2026.1.zip"), overwrite = TRUE)
utils::unzip(file.path(mdir, "msigdb.2026.1.zip"), exdir = mdir)
Sys.setenv(R_USER_CACHE_DIR = cache)

rgz <- function(p) read.delim(if (grepl("gz$", p)) gzfile(p) else p, check.names = FALSE, stringsAsFactors = FALSE)
G <- rgz(arg[["groups"]]); U <- readLines(arg[["universe"]]); fit <- rgz(arg[["fit"]])
if (!setequal(U, fit$gene)) stop("the universe is not the Figure 2C tested genes")
LISTS <- c(split(G$gene, G$group),
           list(G3_CORE_UP_IN_NO_GUT_STUDY = G$gene[G$g3_up_in_no_gut_study %in% c(TRUE, "TRUE")],
                ALL_GUT_EC_HIGHER = fit$gene[fit$call == "NATIVE_EC_HIGHER"], ALL_SC_SIC_HIGHER = fit$gene[fit$call == "SC_SIC_HIGHER"]))
if (!all(unlist(LISTS) %in% U)) stop("a listed gene is outside the tested background")

ENTREZ <- suppressMessages(bitr(U, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
ENTREZ <- ENTREZ[!duplicated(ENTREZ$SYMBOL), ]
HALL <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H"))[, c("gs_name", "gene_symbol")]
P <- list(pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
run <- function(genes) {
  e <- ENTREZ$ENTREZID[ENTREZ$SYMBOL %in% genes]
  r <- list(
    GO_BP = do.call(enrichGO, c(list(gene = genes, universe = U, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP"), P)),
    REACTOME = do.call(enrichPathway, c(list(gene = e, universe = ENTREZ$ENTREZID, organism = "human", readable = TRUE), P)),
    HALLMARK = do.call(enricher, c(list(gene = genes, universe = U, TERM2GENE = HALL), P)))
  do.call(rbind, lapply(names(r), function(k) { d <- if (is.null(r[[k]])) NULL else as.data.frame(r[[k]])
    if (is.null(d) || !nrow(d)) return(NULL)
    data.frame(collection = k, term_id = d$ID, term = d$Description, genes_in_list_and_term = d$Count,
               gene_ratio = d$GeneRatio, bg_ratio = d$BgRatio, p = d$pvalue, p_adj = d$p.adjust, genes = d$geneID) }))
}
RES <- do.call(rbind, lapply(names(LISTS), function(n) {
  r <- run(LISTS[[n]]); if (is.null(r)) NULL else cbind(list = n, list_genes = length(LISTS[[n]]), r) }))
RES$called <- RES$p_adj < 0.05
RES <- RES[order(RES$list, RES$collection, RES$p), ]
D <- RES[RES$called, ]

write_lines <- function(d, con) writeLines(c(paste(names(d), collapse = "\t"), do.call(paste, c(lapply(d, as.character), sep = "\t"))),
                                           con, sep = "\n", useBytes = TRUE)
h <- gzfile(arg[["all_terms"]], "wb"); write_lines(RES, h); close(h)
h <- file(arg[["called"]], "wb"); write_lines(D, h); close(h)
cat(sprintf("pathway enrichment: %d terms tested, %d called\n", nrow(RES), nrow(D)))
