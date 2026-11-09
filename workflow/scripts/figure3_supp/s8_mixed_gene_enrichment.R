# Figure S13A: pathway enrichment of the Figure 3 gene groups with and without their "mixed" genes (genes shared by 2 of
# the 3 intestinal studies). Lists: G1, G2, G3 alone; the mixed genes alone (61 G1-like, 35 G2-like, and the 16 G2-like genes
# higher in SC-SIC that would join G3 if 2 of 3 counted as not shared); and each group with the mixed genes it would join
# (G1 + G1-like, G2 + G2-like, G3 + the 16). Test: one-sided hypergeometric (clusterProfiler) against GO biological process
# (org.Hs.eg.db), Reactome (ReactomePA) and MSigDB Hallmark (MSigDB 2026.1 via msigdbr, read from the release zip);
# Benjamini-Hochberg within list x collection; gene sets of 10-500 genes; called at adjusted p < 0.05. Background: the
# 14,974 genes the donor-level Figure 3D fit tested. Then each group is compared with itself plus its mixed genes:
# every term called in either list, KEPT / LOST / NEW, and the mixed genes inside the term.
# Note: clusterProfiler's order of genes inside a term's gene string ("genes", "mixed_genes_in_term") is not reproducible
# between sessions with the same package versions; the gene sets are.
# Inputs:  groups       gene_groups_definitions.tsv.gz (gene, group)
#          mixed        mixed_genes.tsv (would_be, would_be_g3_if_mixed_counted_not_shared)
#          universe     universe_paneld_tested.txt (background genes)
#          msigdb_zip   msigdb.2026.1.zip (the msigdbr 26.1 release archive)
# Outputs: all          enrichment_all_terms.tsv.gz (every tested term of every list)
#          called       enrichment_called.tsv (terms at adjusted p < 0.05)
#          compare      group_vs_group_with_mixed.tsv (each group against itself with its mixed genes)
#          lists        lists.tsv (list, gene)
#   Rscript s8_mixed_gene_enrichment.R groups=... mixed=... universe=... msigdb_zip=... all=... called=... compare=... lists=...
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA); library(msigdbr) })

# msigdbr 26.1 reads its gene sets from the release zip unpacked in tools::R_user_dir("msigdbr", "cache") (it looks for
# the zip there and the .rds files it holds); set that up in a temporary folder from the given zip, so nothing is downloaded.
msigdb_sets <- function(zip, collection) {
  root <- tempfile("msigdbr_cache_"); Sys.setenv(R_USER_CACHE_DIR = root)
  cache <- tools::R_user_dir("msigdbr", "cache"); dir.create(cache, recursive = TRUE)
  stopifnot(unname(tools::md5sum(zip)) == "512ba99c6827141a9d471972b812d4ac")   # the release MD5 built into msigdbr 26.1
  file.copy(zip, file.path(cache, "msigdb.2026.1.zip")); utils::unzip(zip, exdir = cache)
  as.data.frame(msigdbr(species = "Homo sapiens", collection = collection))[, c("gs_name", "gene_symbol")]
}
rgz <- function(p) read.delim(if (grepl("gz$", p)) gzfile(p) else p, check.names = FALSE, stringsAsFactors = FALSE)

G <- rgz(A[["groups"]]); M <- rgz(A[["mixed"]]); U <- readLines(A[["universe"]])
g3if <- M$would_be_g3_if_mixed_counted_not_shared %in% c(TRUE, "TRUE")
stopifnot(all(M$would_be[g3if] == "G2"))
if (!identical(c(sum(M$would_be == "G1"), sum(M$would_be == "G2"), sum(g3if)), c(61L, 35L, 16L)))
  stop("mixed list sizes are not the expected 61 / 35 / 16")
grp <- function(k) G$gene[G$group == k]
G1 <- grp("G1_EC_GENES_SC_SIC_LACK"); G2 <- grp("G2_EC_GENES_SC_SIC_HAVE"); G3 <- grp("G3_SC_SIC_OWN_PROGRAM")
M1 <- M$gene[M$would_be == "G1"]; M2 <- M$gene[M$would_be == "G2"]; M3 <- M$gene[g3if]
LISTS <- list(G1 = G1, G2 = G2, G3 = G3,
              MIXED_G1_LIKE = M1, MIXED_G2_LIKE = M2, MIXED_G3_IF_NOT_SHARED = M3,
              G1_PLUS_MIXED_G1_LIKE = c(G1, M1), G2_PLUS_MIXED_G2_LIKE = c(G2, M2), G3_PLUS_MIXED_G3_IF_NOT_SHARED = c(G3, M3))
if (any(vapply(LISTS, anyDuplicated, 0L) > 0)) stop("duplicate gene in a list")
if (!all(unlist(LISTS) %in% U)) stop("a listed gene is outside the tested background")

ENTREZ <- suppressMessages(bitr(U, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
ENTREZ <- ENTREZ[!duplicated(ENTREZ$SYMBOL), ]
HALL <- msigdb_sets(A[["msigdb_zip"]], "H")
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

# each group against itself with its mixed genes added
CMP <- do.call(rbind, lapply(list(c("G1", "G1_PLUS_MIXED_G1_LIKE"), c("G2", "G2_PLUS_MIXED_G2_LIKE"), c("G3", "G3_PLUS_MIXED_G3_IF_NOT_SHARED")), function(pr) {
  a <- RES[RES$list == pr[1], ]; b <- RES[RES$list == pr[2], ]
  ka <- paste(a$collection, a$term_id); kb <- paste(b$collection, b$term_id)
  keys <- union(ka[a$called], kb[b$called]); if (!length(keys)) return(NULL)
  ia <- match(keys, ka); ib <- match(keys, kb)
  data.frame(group = pr[1], with_mixed = pr[2], collection = sub(" .*", "", keys), term_id = sub("^\\S+ ", "", keys),
             term = ifelse(is.na(ib), a$term[ia], b$term[ib]),
             p_adj_group = a$p_adj[ia], p_adj_with_mixed = b$p_adj[ib],
             genes_group = a$genes_in_list_and_term[ia], genes_with_mixed = b$genes_in_list_and_term[ib],
             status = ifelse(!is.na(ia) & a$called[ia] %in% TRUE,
                             ifelse(!is.na(ib) & b$called[ib] %in% TRUE, "KEPT", "LOST"), "NEW"),
             mixed_genes_in_term = ifelse(is.na(ib), "", vapply(strsplit(b$genes[ib], "/"), function(x)
               paste(intersect(x, setdiff(LISTS[[pr[2]]], LISTS[[pr[1]]])), collapse = "/"), ""))) }))
CMP <- CMP[order(CMP$group, CMP$status, CMP$p_adj_with_mixed), ]

wtsv <- function(D, p, gz = FALSE) { h <- if (gz) gzfile(p, "wb") else file(p, "wb")
  writeLines(c(paste(names(D), collapse = "\t"), do.call(paste, c(lapply(D, as.character), sep = "\t"))), h, sep = "\n", useBytes = TRUE); close(h) }
wtsv(RES, A[["all"]], gz = TRUE)
wtsv(RES[RES$called, ], A[["called"]])
wtsv(CMP, A[["compare"]])
wtsv(data.frame(list = rep(names(LISTS), lengths(LISTS)), gene = unlist(LISTS, use.names = FALSE)), A[["lists"]])
