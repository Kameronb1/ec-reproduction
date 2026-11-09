# Figure S18 A, S18 B (Singh 2024 column), S19 A and S19 B: the value tables behind the identity-factor supplement.
#
# What it computes, from the selection and per-regulon values (identity_selection.py) and the target-placement tests
# (identity_target_placement.py):
#   gut_study_mean / gut_elmentaite_part  for every regulon of every selected factor, EC minus other enteroendocrine in
#       each intestinal study its network was not trained on (Elmentaite 2020 + 2021: mean of the counted parts, and the parts);
#       note = color class (pass_ec, pass_beta: replicated, by direction; fail: not replicated)
#   targets_G2_of_G2G3  factors with targets higher in SC-SIC: pooled targets in G2 + G3 that are G2, of all in G2 + G3, the count expected
#       from group size, and the lean call ("no lean" or "leans G2/G3, q ...")
#   targets_beta_side   factors with targets lower in SC-SIC: pooled targets among beta-side genes, expected count, q
#   singh_effect / singh_passing_same_direction  Singh 2024 induced EC minus other induced enteroendocrine per scored
#       regulon, and per factor the replicated regulons moving with the factor's SC direction, of those scored
# The Singh rows are also written on their own (singh_out). Checks that the agreeing intestinal-study counts and the G2/G3
# split equal the selection's own records.
# Inputs (name=path): values, parts, factors, tests, fdisp      Outputs: out, singh_out
#   Rscript identity_supplement_values.R values=... parts=... factors=... tests=... fdisp=... out=... singh_out=...
args <- commandArgs(trailingOnly = TRUE); a <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
if (.Platform$OS.type == "windows") stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "English_United States.utf8")))
wtsv <- function(x, p) write.table(x, p, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
rd <- function(k) read.delim(a[[k]], stringsAsFactors = FALSE, check.names = FALSE)
tb <- function(x) as.character(x) %in% c("TRUE", "True", "true")
V <- rd("values"); GP <- rd("parts"); FA <- rd("factors"); FA <- FA[tb(FA$selected), ]; TS <- rd("tests"); FD <- rd("fdisp")
V$pass <- tb(V$program_pass)
stopifnot(setequal(FA$factor, unique(V$factor)), setequal(FA$factor, FD$factor))

# rows: two halves by SC direction, then factor symbol
HALF <- c("Higher in SC-SIC", "Lower in SC-SIC")
FA$half <- factor(ifelse(FA$direction == "SC-SIC higher", HALF[1], HALF[2]), HALF); FA$tf <- sub("\\(.*$", "", FA$factor)
FA <- FA[order(FA$half, FA$tf, FA$factor), ]
FA$y <- as.numeric(factor(FA$factor, rev(FA$factor)))
rowy <- function(f) FA$y[match(f, FA$factor)]
colof <- function(pass, dir) ifelse(!pass, "fail", ifelse(dir == "SC-SIC higher", "pass_ec", "pass_beta"))
spread <- function(d) { d <- d[order(d$factor, d$network), ]
  d$y <- rowy(d$factor) + ave(seq_len(nrow(d)), d$factor, FUN = function(i) if (length(i) == 1) 0 else seq(-.30, .30, length.out = length(i))); d }

# ---- intestinal studies ----------------------------------------------------------------------------------------------
G <- GP[GP$status == "counted", ]
GM <- unique(G[, c("program", "factor", "network", "study", "study_mean")])
stopifnot(anyDuplicated(GM[, c("program", "study")]) == 0)
GM$pass <- V$pass[match(GM$program, V$program)]; GM$dir <- V$direction[match(GM$program, V$program)]
for (s in unique(GM$study)) { k <- GM$study == s
  stopifnot(all(abs(GM$study_mean[k] - V[[paste0("gut_mean__", s)]][match(GM$program[k], V$program)]) <= 1e-12)) }
V$pd <- c(`SC-SIC higher` = 1, `SC-SIC lower` = -1, `NA` = 0)[V$sc_program_direction]   # the regulon's own SC direction
agree <- tapply(sign(GM$study_mean) == V$pd[match(GM$program, V$program)], factor(GM$program, V$program), sum); agree[is.na(agree)] <- 0
stopifnot(all(agree == V$gut_agree))
GM$col <- colof(GM$pass, GM$dir)
GM <- do.call(rbind, lapply(split(GM, GM$study), spread))
EP <- G[G$study == "ELMENTAITE", ]

# ---- G2/G3 split (targets higher in SC-SIC) and beta-side overlap (targets lower in SC-SIC) ---------------------------
tt <- function(u, t) { x <- TS[TS$level == "factor" & TS$unit == u & TS$test == t, ]; stopifnot(nrow(x) == 1); x }
qtxt <- function(q) ifelse(q < 0.001, "q < 0.001", ifelse(q >= 0.995, "q = 1", sprintf("q = %s", signif(q, 2))))
up <- FA[FA$half == HALF[1], ]; lo <- FA[FA$half == HALF[2], ]
BU <- do.call(rbind, lapply(up$factor, function(f) { g2 <- tt(f, "lean_G2"); g3 <- tt(f, "lean_G3"); fd <- FD[FD$factor == f, ]
  stopifnot(g2$targets_in_G2_G3 == g3$targets_in_G2_G3, g2$overlap + g3$overlap == g2$targets_in_G2_G3, fd$targets_in_G2_G3 == g2$targets_in_G2_G3)
  lean <- fd$g2_g3_lean; to3 <- grepl("G3", lean)
  data.frame(factor = f, g2 = g2$overlap, g3 = g3$overlap, n = g2$targets_in_G2_G3, e2 = g2$plain_expected,
             call = ifelse(lean == "no lean", "no lean", sprintf("leans %s, %s", ifelse(to3, "G3", "G2"), qtxt(if (to3) g3$q else g2$q)))) }))
BL <- do.call(rbind, lapply(lo$factor, function(f) { b <- tt(f, "beta_side"); fd <- FD[FD$factor == f, ]
  stopifnot(as.integer(fd$beta_side_overlap) == b$overlap)
  data.frame(factor = f, k = b$overlap, e = b$plain_expected, q = b$q) }))

# ---- Singh 2024 ------------------------------------------------------------------------------------------------------
SG <- V[tb(V$singh_scored), ]; SG$col <- colof(SG$pass, SG$direction); SG <- spread(SG)
same <- sign(SG$singh_effect) == SG$d
SC_ <- data.frame(factor = FA$factor, n = as.integer(tapply(SG$pass, factor(SG$factor, FA$factor), sum)),
                  k = as.integer(tapply(SG$pass & same, factor(SG$factor, FA$factor), sum)))

singh <- rbind(data.frame(part = "singh_effect", factor = SG$factor, program = SG$program, item = "GSE238274__primary", value = SG$singh_effect, count = NA, total = NA, expected = NA, q = NA, note = SG$col),
               data.frame(part = "singh_passing_same_direction", factor = SC_$factor, program = NA, item = "", value = NA, count = SC_$k, total = SC_$n, expected = NA, q = NA, note = ""))
wtsv(rbind(data.frame(part = "gut_study_mean", factor = GM$factor, program = GM$program, item = GM$study, value = GM$study_mean, count = NA, total = NA, expected = NA, q = NA, note = GM$col),
           data.frame(part = "gut_elmentaite_part", factor = EP$factor, program = EP$program, item = EP$part, value = EP$part_effect, count = NA, total = NA, expected = NA, q = NA, note = ""),
           data.frame(part = "targets_G2_of_G2G3", factor = BU$factor, program = NA, item = "G2", value = NA, count = BU$g2, total = BU$n, expected = BU$e2, q = NA, note = BU$call),
           data.frame(part = "targets_beta_side", factor = BL$factor, program = NA, item = "beta_side", value = NA, count = BL$k, total = NA, expected = BL$e, q = BL$q, note = ""),
           singh), a[["out"]])
wtsv(singh, a[["singh_out"]])
cat("rows:", nrow(GM), "gut means,", nrow(EP), "Elmentaite parts,", nrow(BU) + nrow(BL), "target tests,", nrow(SG), "Singh regulons\n")
