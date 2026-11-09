# Figure 3H, step 1: the fetal-versus-adult intestinal EC gene list.
#
# Within each eligible intestinal donor (adult or fetal, at least 10 EC and 10 other EEC cells by primary study
# label), EC is compared with other EEC; the gene list is the EC x adult interaction of a donor-paired
# limma-voom model (~ 0 + donor + ec + ec:adult; filterByExpr, TMM, voom, treat at 1.5-fold, BH < 0.05).
# A gene is kept only if it is also fitted with the same sign when each intestinal study is left out in turn
# (Elmentaite 2020 + 2021 pair, Hickey 2023, Egozi 2021), and the 18 genes that define the cell labels
# are removed (8 adult-high and 11 fetal-high genes).
# Detection floor (the drawn list): a kept gene must also be detected (>= 1 read) in at least 10% of the
# EC cells of the eligible donors of the age where it is higher (adult-high genes: adult donors' EC cells,
# fetal-high genes: fetal donors' EC cells; cells pooled over donors). This leaves 7 adult-high and 9
# fetal-high genes, the list panel H draws (Figure 4F uses the same list). The lists without the floor are
# written too: age_scores.py scores every row with them as a control.
#
# Inputs:  gut_sums= gut_detected= gut_groups=  primary intestinal EEC pseudobulk by unit and primary study label
#          excluded=                            the 18 label-defining marker genes (column gene)
# Outputs: adult= fetal=                        the drawn lists (with the detection floor)
#          adult_no_floor= fetal_no_floor=      the same selection without the detection floor
#          calls=                               every gene significant in the full fit, stable or not, with its EC detection share
#
#   Rscript age_genes.R gut_sums=... gut_detected=... gut_groups=... excluded=... adult=... fetal=...
#           adult_no_floor=... fetal_no_floor=... calls=...

args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
suppressPackageStartupMessages({ library(limma); library(edgeR); library(data.table) })

FLOOR_CELLS <- 10; FDR <- 0.05; LFC_FLOOR <- log2(1.5); DETECTION_FLOOR <- 0.10
E20 <- "SF_ELMENTAITE2020_DEVGUT"; E21 <- "SF_ELMENTAITE2021_GUT_ATLAS"; EGOZI <- "SF_EGOZI2021_DEVELOPMENTAL_EEC"; HICKEY <- "SF_HICKEY2023_ADULT_EEC"
STUDY3 <- list(Elmentaite_pair = c(E20, E21), Hickey_2023 = HICKEY, Egozi_2021 = EGOZI)
age_of <- function(fam, unit) fcase(
  fam == HICKEY, "adult",
  fam == EGOZI, fifelse(grepl("_(21|23)w$", unit), "fetal", "neither"),
  fam == E21, fifelse(startsWith(unit, "F"), "fetal", fifelse(startsWith(unit, "A"), "adult", NA_character_)),
  fam == E20, fifelse(grepl("^[0-9]", unit), "fetal", fifelse(startsWith(unit, "T"), "neither", NA_character_)),
  default = NA_character_)
EXCL18 <- fread(arg[["excluded"]])$gene

s <- read.delim(gzfile(arg[["gut_sums"]]), check.names = FALSE)
G <- fread(arg[["gut_groups"]], colClasses = list(character = "unit"))
stopifnot(identical(names(s)[-1], G$group), !anyDuplicated(s$gene))
M <- as.matrix(s[, -1]); rownames(M) <- s$gene
G[, age := age_of(family_id, unit)]
if (anyNA(G$age)) stop("an intestinal donor has no age group")

eligible <- function(g) {
  w <- dcast(g[, .(n = sum(n_cells)), by = .(family_id, unit, age, class)], family_id + unit + age ~ class, value.var = "n", fill = 0)
  w[age %in% c("adult", "fetal") & pmin(POSITIVE, NEGATIVE) >= FLOOR_CELLS]
}
fit_age <- function(donors) {
  donors <- donors[order(family_id, unit)]
  cols <- c(rbind(paste(donors$family_id, donors$unit, "POSITIVE", sep = "|"), paste(donors$family_id, donors$unit, "NEGATIVE", sep = "|")))
  donor <- factor(rep(paste(donors$family_id, donors$unit, sep = "|"), each = 2))
  ec <- rep(c(1, 0), nrow(donors)); adult <- rep(as.numeric(donors$age == "adult"), each = 2)
  design <- model.matrix(~ 0 + donor + ec + ec:adult); colnames(design)[colnames(design) == "ec:adult"] <- "ec_adult"
  y <- DGEList(M[, cols, drop = FALSE])
  y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  fit <- lmFit(voom(y, design), design)
  tr <- topTreat(treat(fit, lfc = LFC_FLOOR), coef = "ec_adult", number = Inf, sort.by = "none")
  data.table(gene = rownames(y), logFC = tr$logFC, adjP = tr$adj.P.Val,
             call = fifelse(tr$adj.P.Val <= FDR, fifelse(tr$logFC > 0, "ADULT_HIGH", "FETAL_HIGH"), "NS"))
}

D <- eligible(G)
full <- fit_age(D)
ST <- full[call != "NS"]
for (k in names(STUDY3)) {
  r <- fit_age(D[!family_id %in% STUDY3[[k]]]); i <- match(ST$gene, r$gene)
  ST[, paste0("sign_without_", k) := fifelse(is.na(i), "NOT_TESTED", fifelse(sign(r$logFC[i]) == sign(logFC), "same", "flipped"))]
}
ST[, excluded18 := gene %in% EXCL18]
ST[, stable := Reduce(`&`, lapply(paste0("sign_without_", names(STUDY3)), function(k) get(k) == "same"))]

# detection share of each significant gene among the eligible donors' EC cells (pooled) at the age where it is higher
Dt <- read.delim(gzfile(arg[["gut_detected"]]), check.names = FALSE)
stopifnot(identical(names(Dt)[-1], G$group), identical(Dt$gene, s$gene))
DM <- as.matrix(Dt[, -1]); rownames(DM) <- Dt$gene
ecg <- G[class == "POSITIVE" & paste(family_id, unit) %in% paste(D$family_id, D$unit)]
share <- function(age_) { g <- ecg[age == age_]; rowSums(DM[, g$group, drop = FALSE]) / sum(g$n_cells) }
SA <- share("adult"); SF <- share("fetal")
ST[, ec_detected_share := fifelse(call == "ADULT_HIGH", SA[gene], SF[gene])]
ST[, detected_10pct := ec_detected_share >= DETECTION_FLOOR]

lines <- function(x, path) { out <- file(path, "wb"); writeLines(x, out); close(out) }
lines(ST[stable & !excluded18 & detected_10pct & call == "ADULT_HIGH", gene], arg[["adult"]])
lines(ST[stable & !excluded18 & detected_10pct & call == "FETAL_HIGH", gene], arg[["fetal"]])
lines(ST[stable & !excluded18 & call == "ADULT_HIGH", gene], arg[["adult_no_floor"]])
lines(ST[stable & !excluded18 & call == "FETAL_HIGH", gene], arg[["fetal_no_floor"]])
fwrite(ST[order(!stable, call, -abs(logFC))], arg[["calls"]], sep = "\t")
n <- function(f, k) sum(ST$stable & !ST$excluded18 & f & ST$call == k)
cat(sprintf("age fit: %d donors; stable adult-high %d, fetal-high %d; with the detection floor %d and %d\n", nrow(D),
            n(TRUE, "ADULT_HIGH"), n(TRUE, "FETAL_HIGH"), n(ST$detected_10pct, "ADULT_HIGH"), n(ST$detected_10pct, "FETAL_HIGH")))
