# Export the exact MSigDB C8 (cell type signature) members used by the C8 enrichment
# as plain text, with per-set source citations.
#
# Input: the msigdbr 26.1.1 cache built from Zenodo record 18968178
#   (msigdb.2026.1.zip, MD5 512ba99c6827141a9d471972b812d4ac), laid out as
#   <cache_root>/R/msigdbr/{msigdb.2026.1.zip, msigdb.2026.1.Hs.C8.rds, msigdb.2026.1.summary.rds}
# Usage: Rscript export_c8_members.R <cache_root> <out_dir>
#
# The enrichment consumes msigdbr(species = "Homo sapiens", collection = "C8")[, c("gs_name", "gene_symbol")].
# This script writes that exact two-column table (row order preserved) and checks it
# against the raw cached data file, so the plain-text list is the analyzed input, not an approximation.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
cache_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- args[[2]]
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) stop("output directory not empty: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({library(msigdbr); library(data.table); library(digest)})
stopifnot(as.character(packageVersion("msigdbr")) == "26.1.1")

Sys.setenv(R_USER_CACHE_DIR = cache_root)
cache <- normalizePath(tools::R_user_dir("msigdbr", "cache"), winslash = "/", mustWork = TRUE)
expected <- c(
  "msigdb.2026.1.zip" = "84c8ff3270db946ced2a7a4beacf866a3d4ce138af2cbdb1dae76eacb4d4c576",
  "msigdb.2026.1.Hs.C8.rds" = "753e8355fe59dd007fc7509d839c4f866d3d5cf92ae2232a47a9ba6e276de0d9",
  "msigdb.2026.1.summary.rds" = "f943c4ab848a29708145d429f535b2ba574b27cc9cf2ffec4bbadf3da23a7f63")
for (f in names(expected)) {
  got <- digest(file.path(cache, f), algo = "sha256", file = TRUE)
  if (!identical(got, unname(expected[[f]]))) stop("sha256 mismatch for ", f, ": ", got)
}
before <- file.info(list.files(cache, full.names = TRUE))[, c("size", "mtime")]

# Exactly what enrich_c8.R consumes.
C8 <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C8"))
consumed <- C8[, c("gs_name", "gene_symbol")]

# Independent route: the raw cached data frame.
raw <- readRDS(file.path(cache, "msigdb.2026.1.Hs.C8.rds"))
stopifnot(all(raw$db_version == "2026.1.Hs"), all(raw$gs_collection == "C8"))
# msigdbr reorders rows by set name; the (set, symbol) pairs must be the same multiset,
# including pairs repeated because several source genes map to one symbol.
pair_c <- sort(paste(consumed$gs_name, consumed$gene_symbol, sep = "\t"), method = "radix")
pair_r <- sort(paste(raw$gs_name, raw$db_gene_symbol, sep = "\t"), method = "radix")
if (!identical(pair_c, pair_r)) stop("msigdbr() (set, symbol) multiset differs from the raw cached C8 table")
if (!identical(C8$db_gene_symbol, C8$gene_symbol)) stop("human symbols were remapped by msigdbr")

after <- file.info(list.files(cache, full.names = TRUE))[, c("size", "mtime")]
if (!identical(before, after)) stop("msigdbr modified the cache")

# Rows in exactly the order and multiplicity the enrichment receives as TERM2GENE.
members <- data.table(gs_id = C8$gs_id, gs_name = C8$gs_name, gene_symbol = C8$gene_symbol,
                      ncbi_gene = C8$ncbi_gene, ensembl_gene = C8$ensembl_gene)
fwrite(members, file.path(out_dir, "c8_members_msigdb_2026.1.Hs.tsv"), sep = "\t", eol = "\n")

sets <- unique(as.data.table(raw)[, .(gs_id, gs_name, gs_collection, gs_collection_name, gs_source_species,
                                      gs_pmid, gs_geoid, gs_exact_source, gs_url, gs_description)])
stopifnot(!anyDuplicated(sets$gs_id), !anyDuplicated(sets$gs_name))
n <- members[, .(n_rows = .N, n_unique_symbols = uniqueN(gene_symbol)), by = gs_id]
sets <- merge(sets, n, by = "gs_id", sort = FALSE)
fwrite(sets, file.path(out_dir, "c8_set_citations_msigdb_2026.1.Hs.tsv"), sep = "\t", eol = "\n")

writeLines(c(
  paste0("rows\t", nrow(members)),
  paste0("sets\t", nrow(sets)),
  paste0("unique_symbols\t", uniqueN(members$gene_symbol)),
  paste0("repeated_set_symbol_rows\t", sum(duplicated(members[, .(gs_name, gene_symbol)]))),
  paste0("sets_without_pmid\t", sum(sets$gs_pmid == "")),
  paste0("R\t", R.version.string),
  paste0("msigdbr\t", packageVersion("msigdbr")),
  paste0("data.table\t", packageVersion("data.table")),
  paste0("babelgene\t", packageVersion("babelgene"))),
  file.path(out_dir, "export_summary.tsv"))
cat("C8 export complete:", nrow(members), "rows,", nrow(sets), "sets\n")
