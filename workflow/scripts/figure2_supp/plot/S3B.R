# Figure S3B: share of each test study's EC-labeled cells (colored by source) and other terminal EEC cells
# classified as EC as the score cutoff moves, with the cutoff used (dashed), the range of the five leave-one-group-out
# refitted cutoffs carried to this score scale (grey band) and each refitted cutoff (diamonds, EC labels).

CLS <- c("Primary-labeled EC" = "EXPLICIT_EC", "Other terminal EEC" = "TERMINAL_NON_EC_EEC",   # class names as the sweep
         "Unresolved EEC label" = "BROAD_EEC_PARENT_UNRESOLVED")                                     # tables write them
S <- read_table("cutoff_sweep.tsv.gz"); R <- read_table("transported_refits.tsv"); L <- read_table("sweep_landmarks.tsv")
Sh <- read_table("hio_cutoff_sweep.tsv.gz"); Rh <- read_table("hio_transported_refits.tsv")   # McCauley 2023
S$cls <- CLS[S$class]; R$cls <- CLS[R$class]; Sh$cls <- CLS[Sh$class]; Rh$cls <- CLS[Rh$class]
S <- rbind(S, Sh[names(S)]); R <- rbind(R, Rh[names(R)])
live <- unique(L$live_cutoff); band <- c(unique(L$band_lo), unique(L$band_hi))   # cutoff used; range of the refits on this scale
FAMS <- c("gut_cell_atlas", "hickey", "egozi", "gse238274", "gse214852_hio")
SRC_KEY <- c(gut_cell_atlas = "EC_GUT", hickey = "EC_GUT", egozi = "EC_GUT", gse238274 = "EC_INDUCED", gse214852_hio = "EC_HIO")
CLSB <- c(EC_GUT = EC_SRC[["gut"]], EC_INDUCED = EC_SRC[["induced"]], EC_HIO = EC_SRC[["hio"]], TERMINAL_NON_EC_EEC = CLASS[["TERMINAL_NON_EC_EEC"]])
CLSB_LAB <- c(EC_GUT = "EC (native intestinal)", EC_INDUCED = "EC (induced)", EC_HIO = "EC (HIO)", TERMINAL_NON_EC_EEC = "Other terminal EEC")
HEAD_BOLD <- paste0(study(FAMS), ifelse(FAMS == "hickey", " (nuclei)", "")); HEAD_TYPE <- study_type(FAMS)
S2 <- S[S$cls %in% c("EXPLICIT_EC", "TERMINAL_NON_EC_EEC"), ]; S2 <- S2[order(S2$family, S2$cls, S2$cutoff), ]
S2$st <- factor(S2$family, levels = FAMS); S2$key <- ifelse(S2$cls == "EXPLICIT_EC", SRC_KEY[S2$family], S2$cls)
R2 <- R[R$cls == "EXPLICIT_EC", ]; R2$st <- factor(R2$family, levels = FAMS)
H_BUDGET_B <- 2.65; KEY_GAP_B <- .06; KEY_H_B <- .34                 # panel height to fill, gap above the key, key (in)
pb <- ggplot(S2, aes(cutoff, post_rate)) +
  annotate("rect", xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf, fill = "#EDEDED") +
  geom_vline(xintercept = live, linewidth = .35, linetype = "dashed", color = "black") +
  geom_step(aes(color = key, group = cls), linewidth = .5, direction = "hv") +
  geom_point(data = R2, aes(equiv_cutoff, post_rate), shape = 23, fill = "white", color = "black", size = 1.3, stroke = .35) +
  scale_color_manual(values = CLSB, guide = "none") +
  scale_x_continuous(breaks = c(-.5, 0, .5, 1, 1.5), labels = c("-0.5", "0", "0.5", "1", "1.5")) +
  scale_y_continuous(breaks = c(0, 50, 100), labels = function(x) paste0(x, "%")) +
  coord_cartesian(xlim = c(-0.9, 1.5), ylim = c(-2.5, 103), expand = FALSE) +
  facet_grid(. ~ st) +
  labs(x = "Four-marker score cutoff", y = "Called EC (% of cells\nwith that label)") + theme_kit("y") +
  theme(strip.text = element_blank(), panel.spacing.x = unit(GRID_BC$S, "in"), plot.margin = margin(0, 0, 0, 0),
        axis.title.y = element_text(lineheight = .95, margin = margin(r = 3)))
g <- ggplotGrob(pb)
H_P <- floor(100 * (H_BUDGET_B - TOP - HEAD_H - other_h(g) - KEY_GAP_B - KEY_H_B - PAD_B)) / 100   # plot height (in)
g <- fix_rows(fix_columns(g, GRID_BC$G, GRID_BC$P, GRID_BC$S, GRID_BC$W), H_P)
g <- add_headers(g, lapply(seq_along(FAMS), function(i) list(l = i, r = i, bold = HEAD_BOLD[i], plain = HEAD_TYPE[i])))
key_b <- list(list(list(label = sprintf("Cutoff used (%.4f)", live), glyph = "dashv", col = "black"),
                   list(label = "Range of refitted cutoffs, matched to this scale", glyph = "square", col = "#EDEDED"),
                   list(label = "Refitted cutoff, matched to this scale", glyph = "diamond", col = "black")),
              lapply(names(CLSB), function(k) list(label = CLSB_LAB[[k]], glyph = "line", col = CLSB[[k]])))
H_G <- total_h(g)
H_B <- TOP + H_G + KEY_GAP_B + KEY_H_B + PAD_B

draw_panel <- function()
  place(stack_pieces(list(list(h = H_G, draw = g), list(h = KEY_GAP_B, draw = NULL), list(h = KEY_H_B, draw = draw_key(key_b)),
                          list(h = PAD_B, draw = NULL))),
        120, 2003, 4380, 1340, GRID_BC$W, H_B - (.30 + .17))
