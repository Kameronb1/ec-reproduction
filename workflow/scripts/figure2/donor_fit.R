# Figure 2C: genes that differ between SC-SIC cells and primary intestinal EC cells, donor-level limma-voom.
#
# Every primary intestinal donor (EC cells by published label) and every sequenced SC-islet sample (SC-SIC
# cells) is one observation: raw counts summed over its cells. Genes are those measured in every profile,
# minus mitochondrial and ribosomal-protein genes; the tested genes are those kept by edgeR filterByExpr on
# the ten pooled study profiles. TMM, voom, then duplicateCorrelation with study as a block (Elmentaite 2020
# and 2021 form one block; three intestinal blocks, six SC blocks), each run twice as in the limma guide; design
# ~ arm; eBayes; BH; called at FDR < 0.05. The seven Elmentaite 2020 child Crohn's disease donors are removed.
# Also writes the pooled-profile voom log2 CPM used by the Figure 2D gene table.
#
# Inputs:  gut_sums= gut_groups=  primary intestinal EEC pseudobulk by donor and published label (POSITIVE = EC)
#          roster=                intestinal reference roster (the excluded Elmentaite 2020 donors)
#          sc_dir=                folder holding the eleven SC-SIC per-sample raw pseudobulks (gene, count)
# Outputs: levels= (fitted arm levels and calls of every tested gene) profiles= tested=
#
#   Rscript donor_fit.R gut_sums=... gut_groups=... roster=... sc_dir=... levels=... profiles=... tested=...

args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
suppressPackageStartupMessages({ library(jsonlite); library(edgeR); library(limma) })
# gene order follows sort(); an English UTF-8 collation gives the published gene order
for (loc in c("English_United States.utf8", "en_US.UTF-8", "C.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc)))) break

FAMILIES <- c(SF_ELMENTAITE2021_GUT_ATLAS = "GCA2021", SF_HICKEY2023_ADULT_EEC = "HICKEY2023",
              SF_EGOZI2021_DEVELOPMENTAL_EEC = "EGOZI2021", SF_ELMENTAITE2020_DEVGUT = "ELM2020")
# SC-SIC sequenced samples of the six 10x studies: study, sample, cells, file in sc_dir
SC <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
study     sample                  cells file
GSE151117 GSE151117__SRR11837751    858 GSE151117__GSE151117__SRR11837751__SC_SIC.tsv.gz
GSE167880 GSE167880__SRR13805895    392 GSE167880__GSE167880__SRR13805895__SC_SIC.tsv.gz
GSE167880 GSE167880__SRR13805896    222 GSE167880__GSE167880__SRR13805896__SC_SIC.tsv.gz
GSE167880 GSE167880__SRR13805897    229 GSE167880__GSE167880__SRR13805897__SC_SIC.tsv.gz
GSE185036 hCiPSC-islets_S6D2        401 GSE185036_hCiPSC-islets_S6D2__hCiPSC-islets_S6D2__SC_SIC.tsv.gz
GSE202497 GSE202497__SRR19140222    395 GSE202497__GSE202497__SRR19140222__SC_SIC.tsv.gz
GSE202497 GSE202497__SRR19140223    686 GSE202497__GSE202497__SRR19140223__SC_SIC.tsv.gz
GSE226346 S6_WT                     770 GSE226346_MUC18396__S6_WT__SC_SIC.tsv.gz
GSE199636 GSE199636__sc_islet_1    3828 GSE199636__sc_islet_1.tsv.gz
GSE199636 GSE199636__sc_islet_2    3443 GSE199636__sc_islet_2.tsv.gz
GSE199636 GSE199636__s6_week2      5584 GSE199636__s6_week2.tsv.gz")
MAIN <- c("GSE151117", "GSE167880", "GSE185036", "GSE199636", "GSE202497", "GSE226346")
LABEL_GENES <- c("TPH1", "SLC18A1", "FEV", "LMX1A")

roster <- fromJSON(arg[["roster"]])
groups <- read.delim(arg[["gut_groups"]], stringsAsFactors = FALSE, colClasses = c(unit = "character"))
sums <- read.delim(gzfile(arg[["gut_sums"]]), check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(sums$gene) || !identical(names(sums)[-1], groups$group)) stop("intestinal sums and groups disagree")
drop <- groups$family_id == roster$family_id & groups$unit %in% roster$excluded_units
if (!setequal(groups$unit[drop], roster$excluded_units)) stop("a roster unit is missing from the intestinal pseudobulk")
ec <- groups$class == "POSITIVE" & !drop
gut <- as.matrix(sums[, -1]); storage.mode(gut) <- "double"; rownames(gut) <- sums$gene

sc <- lapply(seq_len(nrow(SC)), function(i) {
  a <- read.delim(gzfile(file.path(arg[["sc_dir"]], SC$file[i])), stringsAsFactors = FALSE)
  if (anyDuplicated(a$gene)) stop("duplicate gene in ", SC$file[i])
  setNames(as.numeric(a$count), a$gene)
})
common <- sort(Reduce(intersect, c(list(sums$gene), lapply(sc, names))))
common <- common[!grepl("^MT-|^RP[SL][0-9]", common, ignore.case = TRUE)]
if (!all(LABEL_GENES %in% common)) stop("a label gene is outside the common gene set")

# one column per intestinal donor (ordered by study then donor) and per SC sample (ordered by study then sample)
cols <- list(); units <- list()
for (fam in names(FAMILIES)) {
  j <- which(ec & groups$family_id == fam); j <- j[order(groups$unit[j])]
  for (k in j) { cols[[length(cols) + 1]] <- gut[common, k]
    units[[length(units) + 1]] <- data.frame(arm = "PRIMARY_EC", study = FAMILIES[[fam]], unit = groups$unit[k], cells = groups$n_cells[k]) }
}
for (st in MAIN) for (smp in sort(unique(SC$sample[SC$study == st]))) {
  i <- which(SC$study == st & SC$sample == smp)
  cols[[length(cols) + 1]] <- Reduce(`+`, lapply(sc[i], function(v) v[common]))
  units[[length(units) + 1]] <- data.frame(arm = "SC_SIC", study = st, unit = smp, cells = sum(SC$cells[i]))
}
U <- do.call(rbind, units)
Y <- do.call(cbind, cols); rownames(Y) <- common; colnames(Y) <- paste(U$study, U$unit, sep = "|")
if (anyNA(Y)) stop("a profile lacks a common gene")
if (ncol(Y) != 44L) stop("expected 44 observations (33 intestinal donors, 11 SC samples), found ", ncol(Y))
arm <- factor(U$arm, levels = c("PRIMARY_EC", "SC_SIC"))
blk <- ifelse(U$study %in% c("GCA2021", "ELM2020"), "ELM_PAIR", U$study)

# tested genes: filterByExpr on one pooled profile per study
pooled <- vapply(c(unname(FAMILIES), MAIN), function(s) rowSums(Y[, U$study == s, drop = FALSE]), numeric(nrow(Y)))
rownames(pooled) <- common
parm <- factor(c(rep("PRIMARY_EC", length(FAMILIES)), rep("SC_SIC", length(MAIN))), levels = c("PRIMARY_EC", "SC_SIC"))
pd <- DGEList(counts = pooled, group = parm)
tested <- rownames(pd)[filterByExpr(pd, group = parm)]
pd <- calcNormFactors(pd[tested, , keep.lib.sizes = FALSE], method = "TMM")
pv <- voom(pd, model.matrix(~ parm))
colnames(pv$E) <- c(unname(FAMILIES), MAIN)

# donor-level fit, study as a block
d <- calcNormFactors(DGEList(counts = Y, group = arm)[tested, , keep.lib.sizes = FALSE], method = "TMM")
des <- model.matrix(~ arm)
v1 <- voom(d, des); c1 <- duplicateCorrelation(v1, des, block = blk)$consensus.correlation
v2 <- voom(d, des, block = blk, correlation = c1); c2 <- duplicateCorrelation(v2, des, block = blk)$consensus.correlation
fit <- eBayes(lmFit(v2, des, block = blk, correlation = c2))
tt <- topTable(fit, coef = "armSC_SIC", number = Inf, sort.by = "none")[tested, ]
lev <- data.frame(gene = tested, logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val,
                  B = tt$B, label_gene = tested %in% LABEL_GENES)
lev$call <- ifelse(lev$FDR < 0.05 & lev$logFC > 0, "SC_SIC_HIGHER", ifelse(lev$FDR < 0.05 & lev$logFC < 0, "NATIVE_EC_HIGHER", "NO_CALL"))
co <- fit$coefficients[tested, ]
lev$primary_level <- co[, "(Intercept)"]; lev$sc_level <- co[, "(Intercept)"] + co[, "armSC_SIC"]

wgz <- function(x, path) { h <- gzfile(path, "wt"); on.exit(close(h)); write.table(x, h, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA") }
wgz(lev, arg[["levels"]])
wgz(data.frame(gene = rownames(pv$E), pv$E, check.names = FALSE), arg[["profiles"]])
writeLines(tested, arg[["tested"]])
cat(sprintf("donor fit: %d observations, %d tested genes, within-study correlation %.4f; calls %s\n", ncol(Y), length(tested), c2,
            paste(names(table(lev$call)), table(lev$call), collapse = ", ")))
