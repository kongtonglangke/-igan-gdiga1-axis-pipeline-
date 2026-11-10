## =====================================================================
## 08_figures.R  -  composite MR figure (6 panels) + supplementary figure
## fig/FigureMR_composite.pdf/.png/.tiff ; fig/FigureS_MR_details.pdf/.png
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))
library(scales)

harm <- fread(file.path(DIR_OUT, "02_harmonised.tsv"), sep = "\t", data.table = FALSE)
harm <- harm[harm$ok, ]
est  <- fread(file.path(DIR_OUT, "03_mr_estimates.tsv"), sep = "\t", data.table = FALSE)
dg   <- fread(file.path(DIR_OUT, "04_mr_diagnostics.tsv"), sep = "\t", data.table = FALSE)
colr <- fread(file.path(DIR_OUT, "06_coloc_mqtl_igan.tsv"), sep = "\t", data.table = FALSE)
st   <- fread(file.path(DIR_OUT, "07_steiger.tsv"), sep = "\t", data.table = FALSE)
pw   <- fread(file.path(DIR_OUT, "12b_power_curve.tsv"), sep = "\t", data.table = FALSE)
ss   <- fread(file.path(DIR_OUT, "05_single_snp.tsv"), sep = "\t", data.table = FALSE)

CPG_COL <- c(cg19603390 = "#762A83", cg19473623 = "#BDBDBD", cg17994788 = "#1B7837",
             cg04827551 = "#B2182B", cg16101574 = "#E08214")

## ---------------- Panel A: locus map ----------------
## The five CpGs span 69.3 kb, but two pairs sit only 2.6 kb (1-2) and 7.2 kb (3-4)
## apart: on this axis that is 0.8 and 2.1 mm, less than one triangle. Those pairs
## are split one level up; every site keeps its true GRCh37 x, and the numbers are
## set beside the markers so they can never sit on a symbol or on each other.
cpg_df <- data.frame(pos = CPG_INFO$pos37 / 1e6, lab = CPG_INFO$cpg,
                     idx = as.character(seq_len(nrow(CPG_INFO))))
cpg_df <- cpg_df[order(cpg_df$pos), ]
cpg_df$y  <- c(1.00, 1.46, 1.46, 1.00, 1.00)
cpg_df$dx <- c(-0.010, 0.010, 0.010, -0.010, 0.010)
iv_df  <- data.frame(pos = harm$pos / 1e6, y = 0)
x0 <- min(cpg_df$pos); x1 <- max(cpg_df$pos)
## key lines built from CPG_INFO so the numbering can never drift from the labels
key1 <- paste(sprintf("%d %s", 1:3, CPG_INFO$cpg[1:3]), collapse = "   ")
key2 <- paste(sprintf("%d %s", 4:5, CPG_INFO$cpg[4:5]), collapse = "   ")

pA <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
  geom_segment(aes(x = x0, xend = x1, y = 2.02, yend = 2.02), linewidth = 0.35, colour = "grey30") +
  geom_segment(aes(x = x0, xend = x0, y = 1.95, yend = 2.09), linewidth = 0.35, colour = "grey30") +
  geom_segment(aes(x = x1, xend = x1, y = 1.95, yend = 2.09), linewidth = 0.35, colour = "grey30") +
  annotate("text", x = mean(cpg_df$pos), y = 2.20, label = "69.3 kb", size = 2.2, colour = "grey25") +
  geom_segment(data = cpg_df, aes(x = pos, xend = pos, y = 0.08, yend = y - 0.10),
               linewidth = 0.18, linetype = 3, colour = "grey60") +
  geom_point(data = cpg_df, aes(pos, y), shape = 17, size = 1.05, colour = PAL[["mQTL"]]) +
  geom_point(data = iv_df,  aes(pos, y), shape = 21, size = 1.35,
             fill = PAL[["IgAN"]], colour = "white", stroke = 0.3) +
  geom_text(data = cpg_df, aes(pos + dx, y, label = idx), size = 2.0,
            fontface = "bold", colour = PAL[["mQTL"]]) +
  annotate("text", x = 7.095, y = -0.45, hjust = 0, size = 1.8, colour = "grey20",
           label = key1) +
  annotate("text", x = 7.095, y = -0.85, hjust = 0, size = 1.8, colour = "grey20",
           label = key2) +
  coord_cartesian(xlim = c(7.09, 7.38), ylim = c(-1.05, 2.40), clip = "off") +
  scale_x_continuous(name = "Chromosome 7 position (Mb, GRCh37)", breaks = c(7.1, 7.2, 7.3)) +
  scale_y_continuous(breaks = NULL) +
  labs(title = "C1GALT1 locus",
       subtitle = "5 cis-mQTL CpGs (triangles) and their instruments (circles)") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

## ---------------- Panel B: forest plot ----------------
prim <- est[est$method == "IVW (fixed)" | (est$method == "Wald ratio" & est$k == 1), ]
prim <- prim[!is.na(prim$or), ]
prim$label <- sprintf("%s  (k=%d)", prim$cpg, prim$k)
prim$label <- factor(prim$label, levels = prim$label[order(prim$b)])

pB <- ggplot(prim, aes(or, label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.14, linewidth = 0.45, colour = "grey25") +
  geom_point(size = 2.1, shape = 15, colour = PAL[["mQTL"]]) +
  geom_text(aes(x = 1.44, label = sprintf("%.2f (%.2f-%.2f)", or, lo, hi)),
            size = 1.9, hjust = 0, colour = "grey20") +
  coord_cartesian(xlim = c(0.78, 2.1)) +
  scale_x_log10(name = "OR for IgAN (95% CI), per 1-unit increase in methylation M-value",
                breaks = c(0.8, 1.0, 1.2)) +
  labs(title = "MR estimates", subtitle = "IVW (fixed) or Wald ratio; all P > 0.05") +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 6.2),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

## ---------------- Panel C: instrument-level scatter ----------------
sub <- harm[harm$cpg %in% names(which(table(harm$cpg) >= 2)), ]
sub$facet <- factor(sprintf("%s  (k = %d)", sub$cpg, as.integer(table(harm$cpg)[sub$cpg])))
sl  <- est[est$cpg %in% unique(sub$cpg) & est$method == "Wald ratio", ]
sl$facet <- factor(sprintf("%s  (k = %d)", sl$cpg, sl$k), levels = levels(sub$facet))

pC <- ggplot(sub, aes(beta_exposure, beta_outcome)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_abline(data = sl, aes(slope = b, intercept = 0), linetype = 2,
              linewidth = 0.45, colour = PAL[["mQTL"]]) +
  geom_errorbar(aes(ymin = beta_outcome - 1.96 * se_outcome,
                    ymax = beta_outcome + 1.96 * se_outcome),
                width = 0, linewidth = 0.35, colour = "grey35") +
  geom_errorbarh(aes(xmin = beta_exposure - 1.96 * se_exposure,
                     xmax = beta_exposure + 1.96 * se_exposure),
                 height = 0, linewidth = 0.35, colour = "grey35") +
  geom_point(size = 1.9, colour = PAL[["IgAN"]]) +
  facet_wrap(~ facet, scales = "free", ncol = 1) +
  scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 3) +
  labs(title = "Instrument-level effects",
       subtitle = "SNP-methylation vs SNP-IgAN; dashed line = Wald slope",
       x = expression(beta[mQTL]~"(methylation M-value)"),
       y = expression(beta[IgAN]~"(log OR)"))

## ---------------- Panel D: colocalisation ----------------
cl  <- colr[!is.na(colr$PPH4), c("cpg", "PPH0", "PPH1", "PPH2", "PPH3", "PPH4")]
cll <- pivot_longer(cl, -cpg, names_to = "hyp", values_to = "pp")
cll$hyp <- factor(cll$hyp, levels = c("PPH0", "PPH1", "PPH2", "PPH3", "PPH4"))
cll$cpg <- factor(cll$cpg, levels = CPG_INFO$cpg)

pD <- ggplot(cll, aes(cpg, pp, fill = hyp)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.15) +
  geom_text(data = cl, aes(x = cpg, y = 0.5, label = sprintf("H4\n%.3f", PPH4)),
            inherit.aes = FALSE, colour = "white", size = 1.85, lineheight = 0.95) +
  scale_fill_manual(values = c(PPH0 = "#E5E5E5", PPH1 = "#8DA0CB", PPH2 = "#66C2A5",
                               PPH3 = "#FFD92F", PPH4 = "#E41A1C"), name = NULL, drop = FALSE) +
  scale_y_continuous(name = "Posterior probability", expand = c(0, 0),
                     limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
  scale_x_discrete(labels = function(x) sub("cg", "", x)) +
  labs(title = "Colocalisation, mQTL x IgAN",
       subtitle = sprintf("PPH3 = %.3f-%.3f; PPH4 = %.3f-%.3f",
                          min(cl$PPH3), max(cl$PPH3), min(cl$PPH4), max(cl$PPH4))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
        axis.title.x = element_blank(),
        legend.position = "right", legend.key.size = unit(2.6, "mm"),
        legend.spacing.y = unit(0.4, "mm"))

## ---------------- Panel E: power curves ----------------
pw$cpg <- factor(pw$cpg, levels = CPG_INFO$cpg)
pE <- ggplot(pw, aes(OR, power, colour = cpg)) +
  geom_hline(yintercept = 0.8, linetype = 3, linewidth = 0.3, colour = "grey45") +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = CPG_COL, name = NULL,
                      breaks = c("cg04827551", "cg16101574", "cg17994788", "cg19603390")) +
  scale_y_continuous(name = "Power", limits = c(0, 1.01), breaks = c(0, 0.5, 0.8, 1)) +
  scale_x_continuous(name = "True OR per 1-unit M-value", breaks = seq(1.0, 1.5, 0.1)) +
  labs(title = "Calibrated power", subtitle = "10,000 replicates per cell; dotted line = 80%") +
  theme(legend.position = c(0.99, 0.06), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 5.4), legend.spacing.y = unit(0.2, "mm"))

## ---------------- Panel F: Steiger directionality ----------------
s2 <- st[!is.na(st$R2_exposure), ]
s2l <- rbind(data.frame(cpg = s2$cpg, what = "Methylation", r2 = s2$R2_exposure),
             data.frame(cpg = s2$cpg, what = "IgAN",        r2 = s2$R2_outcome))
s2l$cpg <- factor(s2l$cpg, levels = CPG_INFO$cpg)
s2l$what <- factor(s2l$what, levels = c("Methylation", "IgAN"))

pF <- ggplot(s2l) +
  geom_line(aes(x = r2, y = cpg, group = cpg), colour = "grey65", linewidth = 0.45) +
  geom_point(aes(x = r2, y = cpg, colour = what), size = 2.15) +
  scale_colour_manual(values = c(Methylation = PAL[["mQTL"]], IgAN = PAL[["IgAN"]]), name = NULL) +
  scale_x_log10(name = expression(R^2~"(log scale)"),
                breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
                labels = c("1e-05", "1e-04", "1e-03", "1e-02", "1e-01", "1"),
                limits = c(5e-6, 1.5)) +
  scale_y_discrete(limits = rev(CPG_INFO$cpg[CPG_INFO$cpg %in% s2l$cpg]),
                   labels = function(x) sub("cg", "", x)) +
  labs(title = "Steiger directionality",
       subtitle = "Methylation R2 2-4 orders higher (P < 1e-15)") +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 5.8),
        legend.position = "top", legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"), legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.y = element_blank())

## ---------------- assemble main figure ----------------
design <- "AABB\nCCDD\nEEFF"
fig <- pA + pB + pC + pD + pE + pF +
  plot_layout(design = design, heights = c(1.0, 1.25, 1.0)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(DIR_FIG, "FigureMR_composite.pdf"), fig, width = 180, height = 172,
       units = "mm", device = cairo_pdf)
## the PNG is what gets embedded in the Word file; it is rendered at the same 600 dpi as the
## TIFF so that the two have identical pixel dimensions (the package checker enforces this)
ggsave(file.path(DIR_FIG, "FigureMR_composite.png"), fig, width = 180, height = 172,
       units = "mm", dpi = 600)
ggsave(file.path(DIR_FIG, "FigureMR_composite.tiff"), fig, width = 180, height = 172,
       units = "mm", dpi = 600, compression = "lzw")

## ================= supplementary figure =================
est2 <- est[!is.na(est$or), ]
est2$method <- factor(est2$method,
  levels = rev(c("Wald ratio", "IVW (fixed)", "IVW (random)", "MR-Egger (slope)",
                 "Weighted median", "Penalised median", "Weighted mode", "Maximum likelihood")))
kk <- table(harm$cpg)
estA <- est2[est2$cpg %in% names(kk)[kk >= 2], ]
sA <- ggplot(estA, aes(or, method, colour = cpg)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.24, linewidth = 0.4,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 1.8, shape = 15, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = CPG_COL, name = NULL, breaks = unique(estA$cpg)) +
  scale_x_log10(name = "OR for IgAN (95% CI)") +
  labs(title = "All sensitivity estimators",
       subtitle = "CpGs with >=2 instruments; single-instrument CpGs are in panel b") +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 5.6),
        legend.position = "top", legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 5.4))

ss$cpg <- factor(ss$cpg, levels = CPG_INFO$cpg)
sB <- ggplot(ss, aes(wald_or, reorder(rsid, wald_or), colour = cpg)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = exp(wald_beta - 1.96 * wald_se),
                     xmax = exp(wald_beta + 1.96 * wald_se)),
                 height = 0.15, linewidth = 0.4) +
  geom_point(size = 1.8, shape = 18) +
  scale_colour_manual(values = CPG_COL, name = NULL,
                      breaks = CPG_INFO$cpg[CPG_INFO$cpg %in% ss$cpg],
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_log10(name = "Wald ratio (95% CI)") +
  labs(title = "Single-instrument estimates",
       subtitle = "Each instrument analysed separately") +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 5.8),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.4, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 5.0), legend.spacing.x = unit(1.2, "mm"))

x <- harm[harm$cpg == "cg04827551", ]
loo <- do.call(rbind, lapply(seq_len(nrow(x)), function(j) {
  b <- sum(x$beta_exposure[-j] * x$beta_outcome[-j] / x$se_outcome[-j]^2) /
       sum(x$beta_exposure[-j]^2 / x$se_outcome[-j]^2)
  s <- sqrt(1 / sum(x$beta_exposure[-j]^2 / x$se_outcome[-j]^2))
  data.frame(left_out = x$rsid[j], or = exp(b), lo = exp(b - 1.96 * s), hi = exp(b + 1.96 * s))
}))
b_all <- sum(x$beta_exposure * x$beta_outcome / x$se_outcome^2) /
         sum(x$beta_exposure^2 / x$se_outcome^2)
s_all <- sqrt(1 / sum(x$beta_exposure^2 / x$se_outcome^2))
loo <- rbind(data.frame(left_out = "None (all)", or = exp(b_all),
                        lo = exp(b_all - 1.96 * s_all), hi = exp(b_all + 1.96 * s_all)), loo)
loo$left_out <- factor(loo$left_out, levels = rev(loo$left_out))
sC <- ggplot(loo, aes(or, left_out)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.4, colour = "grey30") +
  geom_point(size = 1.8, shape = 15, colour = PAL[["mQTL"]]) +
  scale_x_log10(name = "OR for IgAN (95% CI)") +
  labs(title = "Leave-one-out, cg04827551",
       subtitle = "No single instrument drives the estimate") +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 6))

fs <- sA + sB + sC + plot_layout(design = "AAB\nCCB", widths = c(1.22, 1)) +
  plot_annotation(tag_levels = "a")
ggsave(file.path(DIR_FIG, "FigureS_MR_details.pdf"), fs, width = 180, height = 150,
       units = "mm", device = cairo_pdf)
ggsave(file.path(DIR_FIG, "FigureS_MR_details.png"), fs, width = 180, height = 150,
       units = "mm", dpi = 400)

cat("[08] figures written:", paste(list.files(DIR_FIG), collapse = ", "), "\n")
