# Figure S3A: the cells of each test study (Elmentaite 2021 fetal and adult, Hickey 2023, Egozi 2021, Singh 2024)
# and of the McCauley 2023 HIO grafts in that study's own graph embedding, the same cells twice: colored by whether
# the four-marker classifier classifies them as EC (top row) and by primary study label class (bottom row, EC labels
# colored by the cells' source). Each map gives its counts of cells classified or labeled EC; keys under each row.

read_coords <- function(name) read.delim(result_path(name), quote = "", check.names = FALSE, fileEncoding = "UTF-8")
X <- read_coords("heldout_space_coords.tsv.gz")
XH <- read_coords("hio_heldout_space_coords.tsv.gz")            # McCauley 2023 HIO grafts
X <- rbind(X, XH[names(X)])
ARMS <- c("E-MTAB-9536", "E-MTAB-9543", "HICKEY_2023_DRYAD", "EGOZI_2021_ZENODO", "GSE238274", "GSE214852_WT")
unit_of_arm <- unlist(tapply(X$unit_id, X$arm_id, unique))[ARMS]
STUDY <- setNames(study(unit_of_arm), ARMS); TYPE <- setNames(study_type(unit_of_arm), ARMS)
cnt <- data.frame(arm_id = ARMS, cells = as.integer(table(X$arm_id)[ARMS]),
                  called = as.integer(tapply(X$presentation_class == "HIGH_CONFIDENCE_EC", X$arm_id, sum)[ARMS]),
                  ec_label = as.integer(tapply(X$author_label_class == "EXPLICIT_EC", X$arm_id, sum)[ARMS]))
X$call <- factor(ifelse(X$presentation_class == "HIGH_CONFIDENCE_EC", "EC call", "No EC call"), levels = c("EC call", "No EC call"))
X$cls <- factor(X$author_label_class, levels = names(CLASS))
CLS3 <- c(EC_GUT = EC_SRC[["gut"]], EC_INDUCED = EC_SRC[["induced"]], EC_HIO = EC_SRC[["hio"]], CLASS[c("TERMINAL_NON_EC_EEC", "BROAD_EEC_PARENT_UNRESOLVED", "ENDOCRINE_PROGENITOR")])
CLS3_LAB <- c(EC_GUT = "EC (native intestinal)", EC_INDUCED = "EC (induced)", EC_HIO = "EC (HIO)", CLASS_LABEL[c("TERMINAL_NON_EC_EEC", "BROAD_EEC_PARENT_UNRESOLVED", "ENDOCRINE_PROGENITOR")])
k3 <- as.character(X$cls); isec <- k3 == "EXPLICIT_EC"                # EC labels colored by source
k3[isec] <- ifelse(X$arm_id[isec] == "GSE238274", "EC_INDUCED", ifelse(X$arm_id[isec] == "GSE214852_WT", "EC_HIO", "EC_GUT"))
X$cls3 <- factor(k3, levels = names(CLS3))
CALL_COL <- c("EC call" = CALLS[["HIGH_CONFIDENCE_EC"]], "No EC call" = "#D9D9D9")
cnt$top <- sprintf("%s of %s\ncalled EC", cnt$called, fmt(cnt$cells))
cnt$bot <- sprintf("%s of %s\nlabeled EC", cnt$ec_label, fmt(cnt$cells))
X$col <- factor(X$arm_id, levels = ARMS); cnt$col <- factor(cnt$arm_id, levels = ARMS)

W_A <- 7.3; H_BUDGET <- 3.55                                          # panel width and height to fill (in)
S_A <- .07; S_BREAK <- .17; R_A <- .03; KEY1_H <- .21; KEY2_H <- .34  # gaps between maps, right margin, key bands (in)
LAB_PT <- 7.5; LAB_LH <- .92
LAB_BAND <- 2 * LAB_PT * LAB_LH / 72 + .06                            # two count lines and a gap above the points
base <- function(d, ylim_top) ggplot(d, aes(disp_1, disp_2)) +
  coord_cartesian(xlim = c(-1.06, 1.06), ylim = c(-1.06, ylim_top), expand = FALSE, clip = "off") +
  facet_grid(. ~ col) + theme_kit("none") +
  scale_size_manual(values = c(`FALSE` = .5, `TRUE` = .3), guide = "none") +   # smaller points for Hickey 2023's many cells
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 90, face = "bold", size = 8, margin = margin(r = 3)),
        strip.text = element_blank(), panel.spacing.x = unit(S_A, "in"), legend.position = "none",
        plot.margin = margin(0, 0, 0, 0))
o1 <- X[order(X$call != "No EC call"), ]                              # cells classified as EC drawn last
o2 <- X[order(c(BROAD_EEC_PARENT_UNRESOLVED = 1, TERMINAL_NON_EC_EEC = 2, ENDOCRINE_PROGENITOR = 3, EXPLICIT_EC = 4)[as.character(X$cls)]), ]
make_rows <- function(ylim_top) {
  lab_y <- ylim_top - .01 * (ylim_top + 1.06)
  top <- base(o1, ylim_top) + geom_point(aes(color = call, size = arm_id == "HICKEY_2023_DRYAD"), stroke = 0) +
    geom_text(data = cnt, aes(x = -1.055, y = lab_y, label = top), hjust = 0, vjust = 1, lineheight = LAB_LH,
              size = LAB_PT / .pt, family = FONT, inherit.aes = FALSE) +
    scale_color_manual(values = CALL_COL, guide = "none") + labs(y = "Four-marker classifier")
  bot <- base(o2, ylim_top) + geom_point(aes(color = cls3, size = arm_id == "HICKEY_2023_DRYAD"), stroke = 0) +
    geom_text(data = cnt, aes(x = -1.055, y = lab_y, label = bot), hjust = 0, vjust = 1, lineheight = LAB_LH,
              size = LAB_PT / .pt, family = FONT, inherit.aes = FALSE) +
    scale_color_manual(values = CLS3, guide = "none") + labs(y = "Primary study label")
  list(top = ggplotGrob(top), bot = ggplotGrob(bot))
}
# a first build measures the gtables' fixed parts; the map size and the room above the points follow from them
g0 <- make_rows(1.5)
left_need <- function(g) { pc <- panel_cols(g); sum(vapply(seq_len(pc[1] - 1), function(i) in_w(g$widths[i]), 0)) }
G_A <- ceiling(100 * max(left_need(g0$top), left_need(g0$bot))) / 100 + .01
GAPS <- c(S_A, S_A, S_A, S_BREAK, S_A)                               # wider gap between the tissue studies and Singh 2024
P_A <- (W_A - G_A - sum(GAPS) - R_A) / 6
H_U <- floor(100 * (H_BUDGET - TOP - HEAD_H - KEY1_H - KEY2_H - PAD_B - other_h(g0$top) - other_h(g0$bot)) / 2) / 100
Y_TOP <- 1.04 + LAB_BAND * 2.10 / (H_U - LAB_BAND)                    # points fill the panel below the count lines
g <- make_rows(Y_TOP)
g$top <- fix_rows(fix_columns(g$top, G_A, P_A, GAPS, W_A), H_U)
g$bot <- fix_rows(fix_columns(g$bot, G_A, P_A, GAPS, W_A), H_U)
heads <- c(list(list(l = 1, r = 2, bold = STUDY[[1]], rule = TRUE,     # one Elmentaite 2021 header over its fetal and adult maps
                     plain = list(list(x = 0, s = TYPE[[1]]), list(x = P_A + GAPS[1], s = TYPE[[2]])))),
           lapply(3:6, function(i) list(l = i, r = i, bold = paste0(STUDY[[i]], if (ARMS[i] == "HICKEY_2023_DRYAD") " (nuclei)" else ""),
                                        plain = TYPE[[i]])))
g$top <- add_headers(g$top, heads)
key1 <- list(list(list(label = "EC call", glyph = "point", col = CALL_COL[["EC call"]]),
                  list(label = "No EC call", glyph = "point", col = CALL_COL[["No EC call"]])))
key2 <- list(lapply(c("EC_GUT", "EC_INDUCED", "EC_HIO"), function(k) list(label = CLS3_LAB[[k]], glyph = "point", col = CLS3[[k]])),
             lapply(c("TERMINAL_NON_EC_EEC", "BROAD_EEC_PARENT_UNRESOLVED", "ENDOCRINE_PROGENITOR"),
                    function(k) list(label = CLS3_LAB[[k]], glyph = "point", col = CLS3[[k]])))
H_TOP <- total_h(g$top); H_BOT <- total_h(g$bot)
H_A <- TOP + H_TOP + KEY1_H + H_BOT + KEY2_H + PAD_B

draw_panel <- function()
  place(stack_pieces(list(list(h = H_TOP, draw = g$top), list(h = KEY1_H, draw = draw_key(key1)),
                          list(h = H_BOT, draw = g$bot), list(h = KEY2_H, draw = draw_key(key2, grid = TRUE)),
                          list(h = PAD_B, draw = NULL))),
        120, 59, 4380, 1878, W_A, H_A - (.30 + .17))
