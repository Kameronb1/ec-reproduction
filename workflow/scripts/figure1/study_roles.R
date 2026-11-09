# Figure 1B: which studies play which role in the analyses of Figures 2-5.
# The input is the curated role table: one row per study with its biological context, donor ages, and for each of
# Figures 2-5 the panels it appears in (f2..f5) and its role there (c2..c5: train, test, split, shown). This script
# makes two adjustments before the table is drawn:
#   * the paper-wide study name "Yu 2021" (the curated table says "Yu Q 2021");
#   * Figure 2 entries are split by role, so that Fig. 2D (shown or compared) is separate from a test:
#       "A,D|B-E" (split)  -> "A|B,C,E|D"  train|test|shown
#       "B-E"     (test)   -> "B,C,E|D"    test|shown
#       "A,D"     (train)  -> "A|D"        train|shown
# Input:  roles   curated role table (16 studies)
# Output: out     role table as drawn
#   Rscript study_roles.R roles=<study_roles.tsv> out=<study_roles_drawn.tsv>
args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
stopifnot(all(c("roles", "out") %in% names(arg)))

tab <- read.delim(arg[["roles"]], colClasses = "character", na.strings = character(0), encoding = "UTF-8")
stopifnot(nrow(tab) == 16)
for (cc in c("study", "label")) tab[[cc]] <- gsub("Yu Q 2021", "Yu 2021", tab[[cc]], fixed = TRUE)

F2 <- list("A,D|B–E" = c("A|B,C,E|D", "train|test|shown"), "B–E" = c("B,C,E|D", "test|shown"), "A,D" = c("A|D", "train|shown"))
stopifnot(all(tab$f2[nzchar(tab$f2)] %in% names(F2)))
for (k in which(nzchar(tab$f2))) { v <- F2[[tab$f2[k]]]; tab$f2[k] <- v[1]; tab$c2[k] <- v[2] }
write.table(tab, arg[["out"]], sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
