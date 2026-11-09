# Figure 4F: the fetal-versus-adult gene list with the 10% detection floor (the list Figure 3H also takes).
#
# Starts from Figure 3H's list without the floor (8 adult-high and 11 fetal-high genes) and keeps a gene only if it is
# detected (>= 1 read) in at least 10% of the EC cells of the eligible intestinal donors of the age where it is higher
# (adult-high genes: the adult donors' EC cells; fetal-high genes: the fetal donors' EC cells; cells pooled over donors).
# Eligible donors are those of Figure 3H's gene fit: adult or fetal, with at least 10 EC and 10 other EEC cells by primary
# study label. Genes stay in list order.
#
# Inputs: adult=, fetal=<Figure 3H lists without the floor>, gut_detected=<cells detecting each gene per donor and label>,
#         gut_groups=<group table: study, donor, EC POSITIVE / other EEC NEGATIVE, cells>
# Outputs: adult_floor=, fetal_floor=<the lists with the floor>
args <- commandArgs(trailingOnly = TRUE)
A <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
FLOOR_CELLS <- 10; DETECTION_FLOOR <- 0.10
E20 <- "SF_ELMENTAITE2020_DEVGUT"; E21 <- "SF_ELMENTAITE2021_GUT_ATLAS"; EGOZI <- "SF_EGOZI2021_DEVELOPMENTAL_EEC"; HICKEY <- "SF_HICKEY2023_ADULT_EEC"
age_of <- function(fam, unit) ifelse(fam == HICKEY, "adult",
  ifelse(fam == EGOZI, ifelse(grepl("_(21|23)w$", unit), "fetal", "neither"),
  ifelse(fam == E21, ifelse(startsWith(unit, "F"), "fetal", ifelse(startsWith(unit, "A"), "adult", NA)),
  ifelse(fam == E20, ifelse(grepl("^[0-9]", unit), "fetal", ifelse(startsWith(unit, "T"), "neither", NA)), NA))))

L <- list(adult = scan(A[["adult"]], what = "", quiet = TRUE), fetal = scan(A[["fetal"]], what = "", quiet = TRUE))
G <- read.delim(A[["gut_groups"]], colClasses = c(unit = "character"), stringsAsFactors = FALSE)
G$age <- age_of(G$family_id, G$unit); stopifnot(!anyNA(G$age))
D <- read.delim(gzfile(A[["gut_detected"]]), check.names = FALSE); stopifnot(identical(names(D)[-1], G$group))
DM <- as.matrix(D[, -1]); rownames(DM) <- D$gene
donor <- paste(G$family_id, G$unit)
n_ec <- tapply(G$n_cells[G$class == "POSITIVE"], donor[G$class == "POSITIVE"], sum)
n_other <- tapply(G$n_cells[G$class == "NEGATIVE"], donor[G$class == "NEGATIVE"], sum)
ok <- names(n_ec)[names(n_ec) %in% names(n_other)]; ok <- ok[pmin(n_ec[ok], n_other[ok]) >= FLOOR_CELLS]
ec <- G$class == "POSITIVE" & donor %in% ok & G$age %in% c("adult", "fetal")
share <- function(age) { k <- ec & G$age == age; rowSums(DM[, G$group[k], drop = FALSE]) / sum(G$n_cells[k]) }
SH <- list(adult = share("adult"), fetal = share("fetal"))
stopifnot(all(unlist(L) %in% rownames(DM)))
KEEP <- lapply(c(adult = "adult", fetal = "fetal"), function(s) L[[s]][SH[[s]][L[[s]]] >= DETECTION_FLOOR])
for (s in names(KEEP)) { con <- file(A[[paste0(s, "_floor")]], "wb"); writeLines(KEEP[[s]], con); close(con) }
cat(sprintf("detection floor: %d eligible donors (%d adult, %d fetal); adult-high %d of %d genes, fetal-high %d of %d\n",
            length(unique(donor[ec])), length(unique(donor[ec & G$age == "adult"])), length(unique(donor[ec & G$age == "fetal"])),
            length(KEEP$adult), length(L$adult), length(KEEP$fetal), length(L$fetal)))
