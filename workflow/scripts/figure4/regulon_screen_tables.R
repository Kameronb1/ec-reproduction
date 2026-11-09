# Figure S16B, C and D: regulons enriched for the shared graft decline (B) and for cholesterol synthesis (C), and the
# transcription factors whose regulons are hits for the decline in two or more networks (D).
#
# Input is the screen of 1,113 regulons (predicted target sets of one transcription factor in one of 12 gene
# regulatory networks) against two gene sets, read from the screen's "down_soupx" family: the genes that fall by
# >= 0.5 log2 in SC-SIC, SC-beta and SC-alpha in >= 6 of the 8 graft comparisons (SoupX-corrected, unthinned), and
# Reactome cholesterol biosynthesis. Hit = ordinary and expression-matched q both < 0.05 and fold > 1.
#   B, C  the 12 hit regulons with the highest matched fold, then (B) every further NFE2L1, SOX4 or SOX11 hit;
#         fold interval = exact Poisson interval on the targets found divided by the matched expectation
#   D     factors whose regulons are decline hits in >= 2 networks, plus NFE2L1; one tile per factor and network from
#         the activated-target (+) regulon only (the hit if there is one, otherwise the highest matched fold);
#         "hit in N of M": M = networks with a (+) regulon for that factor; networks without one are "no list"
#
# Inputs: screen=<screen_all_lists.tsv.gz>
# Outputs: a=<B table>  b=<C table>  d=<D table>
args <- commandArgs(trailingOnly = TRUE)
A_ <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages(library(data.table))

# network (bank) -> study name and type, in the order D draws them
NET <- data.table(bank = c("GSE151117", "GSE199636", "GSE167880", "GSE185036", "GSE226346", "GSE202497", "GSE221156", "GSE148073",
                           "EGOZI_2021_ZENODO", "support_small_adult_intestine", "E-MTAB-9536", "HICKEY_2023_DRYAD"),
                  net = c("Augsornworawat 2020", "Augsornworawat 2023", "Balboa 2022", "Du 2022", "Siehler 2024", "Zhu 2023", "Bandesh 2026",
                          "Fasolino 2022", "Egozi 2021", "Elmentaite 2021 adult", "Elmentaite 2021 fetal", "Hickey 2023"),
                  group = c(rep("SC-islet", 6), rep("adult islet", 2), rep("intestine", 4)))
NETS <- NET$net

DEC <- "down_soupx"; CHO <- "BIO_PATHWAY__12"
S <- as.data.table(read.delim(gzfile(A_[["screen"]]), stringsAsFactors = FALSE))[family == DEC]
stopifnot(setequal(NET$bank, unique(S$bank)), uniqueN(S$regulon) == 1113L)
S[, net := NET$net[match(bank, NET$bank)]]
S[, fold_lo := ifelse(shared_count == 0, 0, qchisq(0.025, 2 * shared_count) / 2) / matched_expected]
S[, fold_hi := (qchisq(0.975, 2 * (shared_count + 1)) / 2) / matched_expected]
S[, screen_hit := ordinary_q < 0.05 & matched_q < 0.05 & fold > 1]
stopifnot(identical(S$screen_hit, S$hit))

FOCUS <- c("NFE2L1", "SOX4", "SOX11"); N_TOP <- 12L
pick <- function(p) {
  x <- S[partner == p & screen_hit == TRUE][order(rank)]
  top <- x[seq_len(min(N_TOP, .N))]; extra <- x[tf %in% FOCUS & !(regulon %in% top$regulon)]
  rbind(top[, part := "top"], extra[, part := "extra"])
}
GRP_LEV <- c("NFE2L1", "SOX4, SOX11 (also fall in culture alone, where scored)", "other regulons")
grp_of <- function(tf) ifelse(tf %in% c("SOX4", "SOX11"), GRP_LEV[2], ifelse(tf == "NFE2L1", GRP_LEV[1], GRP_LEV[3]))
TIMES <- intToUtf8(0xD7); MINUS <- intToUtf8(0x2212)
fmt_q <- function(q) formatC(q, format = "e", digits = 1)
q_drawn <- function(q) { s <- fmt_q(q); sprintf("%s %s 10^%s%d", sub("e.*$", "", s), TIMES, MINUS, -as.integer(sub("^.*e", "", s))) }
forest <- function(p) { x <- pick(p)
  x[, .(row = seq_len(.N), part, regulon, network = net, bank, tf, sign, group = grp_of(tf),
        label = sprintf("%s%s (%s, %d targets)", tf, ifelse(sign == "-", " repressed", ""), net, measured_targets),
        shared_count, measured_targets, matched_fold, fold_lo, fold_hi, ordinary_q, matched_q, matched_q_drawn = q_drawn(matched_q), rank,
        nfe2l1_band = tf == "NFE2L1")] }

# D: factors hit in >= 2 networks, plus NFE2L1; (+) regulons only
R <- S[partner == DEC & screen_hit == TRUE, .(nets = uniqueN(bank)), by = tf]
CTF <- union(R[nets >= 2][order(-nets, tf), tf], "NFE2L1")
CR <- S[partner == DEC & tf %in% CTF & sign == "+"][order(-screen_hit, -matched_fold)][, .SD[1], by = .(tf, bank)]
NH <- CR[, .(nh = sum(screen_hit)), by = tf][match(CTF, tf)]
NM <- CR[, .(m = uniqueN(bank)), by = tf][match(CTF, tf)]
EMPTY <- CJ(tf = CTF, net = NETS)[!CR, on = .(tf, net)]
D <- rbind(CR[, .(factor = tf, network = net, network_group = NET$group[match(bank, NET$bank)], cell = ifelse(screen_hit, "hit", "tested, not a hit"), regulon, sign,
                  matched_fold, shared_count, measured_targets, tile_top = formatC(matched_fold, format = "f", digits = 1),
                  tile_bottom = sprintf("%d/%d", shared_count, measured_targets))],
           EMPTY[, .(factor = tf, network = net, network_group = NET$group[match(net, NET$net)], cell = "no list", regulon = "", sign = "", matched_fold = NA_real_,
                     shared_count = NA_integer_, measured_targets = NA_integer_, tile_top = "", tile_bottom = "")])
D[, factor := factor(factor, CTF)]; D[, network := factor(network, NETS)]; setorder(D, factor, network)
D <- merge(D, data.table(factor = factor(CTF, CTF), hit_networks = NH$nh, networks_with_list = NM$m,
                         factor_label = sprintf("hit in %d of %d", NH$nh, NM$m)), by = "factor")
setorder(D, factor, network)

wtsv <- function(d, path) { d <- as.data.frame(d)
  for (j in names(d)) { v <- d[[j]]; d[[j]] <- if (is.double(v)) ifelse(is.na(v), "", sprintf("%.15g", v)) else ifelse(is.na(v), "", as.character(v)) }
  con <- file(path, "wb"); on.exit(close(con))
  writeLines(enc2utf8(c(paste(names(d), collapse = "\t"), apply(d, 1, paste, collapse = "\t"))), con, sep = "\n", useBytes = TRUE) }
wtsv(forest(DEC), A_[["a"]]); wtsv(forest(CHO), A_[["b"]]); wtsv(D, A_[["d"]])
