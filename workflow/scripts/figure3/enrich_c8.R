# Figure 3G (cell types): over-representation of MSigDB C8 cell-type signature gene sets in the Figure 3F
# gene groups, and each called term's class.
#
# Same test as the pathway half of 3G: one-sided hypergeometric (clusterProfiler::enricher) against the genes
# the Figure 3D fit tested; sets of 10-500 background genes; Benjamini-Hochberg within each list; called at
# adjusted p < 0.05. Lists: G1, G2, G3, and all Figure 3D calls on each side. Each term is then put in one
# class by the rules below (first match wins): intestinal EC or enteroendocrine, neuroendocrine of other organs,
# glia or neural progenitor, neuron, other. C8 from MSigDB 2026.1 read by msigdbr.
#
# Inputs:  groups= universe= fit=  (gene_groups.py, donor_fit.R)
#          msigdb_zip=              msigdb.2026.1.zip (Zenodo record 18968178)
# Outputs: all_terms= called=      every tested term / the called terms, with term_class
#
#   Rscript enrich_c8.R groups=... universe=... fit=... msigdb_zip=... all_terms=... called=...

args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
suppressPackageStartupMessages({ library(data.table); library(clusterProfiler); library(msigdbr) })

cache <- file.path(tempdir(), "msigdbr_cache"); mdir <- file.path(cache, "R", "msigdbr")
dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
if (unname(tools::md5sum(arg[["msigdb_zip"]])) != "512ba99c6827141a9d471972b812d4ac") stop("msigdb.2026.1.zip has the wrong checksum")
file.copy(arg[["msigdb_zip"]], file.path(mdir, "msigdb.2026.1.zip"), overwrite = TRUE)
utils::unzip(file.path(mdir, "msigdb.2026.1.zip"), exdir = mdir)
Sys.setenv(R_USER_CACHE_DIR = cache)

U <- readLines(arg[["universe"]]); U <- U[nzchar(U)]
DEF <- as.data.table(read.delim(gzfile(arg[["groups"]]), stringsAsFactors = FALSE))
LV <- as.data.table(read.delim(gzfile(arg[["fit"]]), stringsAsFactors = FALSE))
LISTS <- c(split(DEF$gene, DEF$group),
           list(ALL_SC_SIC_HIGHER = LV[call == "SC_SIC_HIGHER", gene], ALL_GUT_EC_HIGHER = LV[call == "NATIVE_EC_HIGHER", gene]))
ORDER <- c("G1_EC_GENES_SC_SIC_LACK", "G2_EC_GENES_SC_SIC_HAVE", "G3_SC_SIC_OWN_PROGRAM", "ALL_SC_SIC_HIGHER", "ALL_GUT_EC_HIGHER")
if (!setequal(names(LISTS), ORDER) || !all(unlist(LISTS) %in% U)) stop("unexpected lists or a gene outside the background")

C8 <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C8"))[, c("gs_name", "gene_symbol")]
P <- list(pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
RES <- rbindlist(lapply(ORDER, function(n) {
  r <- do.call(enricher, c(list(gene = LISTS[[n]], universe = U, TERM2GENE = C8), P)); d <- if (is.null(r)) NULL else as.data.frame(r)
  if (is.null(d) || !nrow(d)) return(NULL)
  data.table(list = n, list_genes = length(LISTS[[n]]), collection = "C8_CELL_TYPE", term_id = d$ID, term = d$Description,
             genes_in_list_and_term = d$Count, gene_ratio = d$GeneRatio, bg_ratio = d$BgRatio, p = d$pvalue, p_adj = d$p.adjust, genes = d$geneID) }))
RES[, called := p_adj < 0.05]; setorder(RES, list, p)
# the table is written as text and read back before classing, so its numbers have the precision of the published table
tmp <- tempfile(fileext = ".tsv"); fwrite(RES, tmp, sep = "\t"); A <- fread(tmp)

GUT_EC <- c("BUSSLINGER_DUODENAL_EC_CELLS", "GAO_LARGE_INTESTINE_ADULT_CA_ENTEROENDOCRINE_CELLS", "DESCARTES_FETAL_INTESTINE_CHROMAFFIN_CELLS")
RULES <- list(
  "Gut EC or enteroendocrine" = function(t) t %in% GUT_EC,
  "Neuroendocrine, other organs" = function(t) grepl("NEUROENDOCRINE|_NE_PRECURSOR", t),
  "Glia or neural progenitor" = function(t) grepl("_OPC$|_HOPC$|_HRGL|OLIGODENDROCYTES|MICROGLIA|_HNPROG$", t),
  "Neuron" = function(t) grepl("NEURON|SYMPATHOBLASTS|MANNO_MIDBRAIN_NEUROTYPES_H(DA|DA2|SERT|OMTN|NBML1|NBML5)$|PURKINJE", t),
  "Other" = function(t) rep(TRUE, length(t)))
cls <- function(t) { out <- rep(NA_character_, length(t)); for (k in names(RULES)) { hit <- is.na(out) & RULES[[k]](t); out[hit] <- k }; out }
A[, term_class := factor(cls(term), names(RULES))]
fwrite(A, arg[["all_terms"]], sep = "\t"); fwrite(A[called == TRUE], arg[["called"]], sep = "\t")
cat(sprintf("C8 enrichment: %d terms tested, %d called\n", nrow(A), sum(A$called)))
