# 在项目根目录独立运行；输入为本次analysis目录中对应分析结果；输出为本图及源数据。
source("scripts/figures/common.R")
setup(3)
a <- list(net = w79$network)
hc <- a$net$dendrograms[[1]]
if (inherits(hc, "dendrogram")) hc <- as.hclust(hc)
nn <- length(hc$order)
pos <- numeric(nn)
pos[hc$order] <- seq_len(nn)
xx <- yy <- numeric(nn - 1)
seg <- list()
hang <- .1
max_h <- max(hc$height)
for (i in seq_len(nn - 1)) {
  kids <- hc$merge[i, ]
  kx <- vapply(kids, function(k) if (k < 0) pos[-k] else xx[k], numeric(1))
  ky <- vapply(kids, function(k) if (k < 0) max(0, hc$height[i] - hang * max_h) else yy[k], numeric(1))
  xx[i] <- mean(kx)
  yy[i] <- hc$height[i]
  seg[[i]] <- rbind(data.frame(x = kx, y = ky, xend = kx, yend = yy[i]), data.frame(x = kx[1], y = yy[i], xend = kx[2], yend = yy[i]))
}
dd <- do.call(rbind, seg)
wr(dd, "C_tree_segments")
pt <- ggplot(dd) +
  geom_segment(aes(x, y, xend = xend, yend = yend), linewidth = .04, colour = "#111111", alpha = .88, lineend = "butt") +
  scale_x_continuous(limits = c(.5, nn + .5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .015))) +
  labs(x = NULL, y = "Height") +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.y = element_line(linewidth = .3, colour = "#303030"),
    axis.ticks.y = element_line(linewidth = .3, colour = "#303030"),
    plot.margin = margin(4, 4, 0, 4)
  )
mm <- data.frame(x = seq_len(nn), module = factor(a$net$colors[a$net$blockGenes[[1]]][hc$order]))
wr(mm, "C_module_order")
bb <- ggplot(mm, aes(x, 1, fill = module)) +
  geom_tile(width = 1, height = 1) +
  scale_x_continuous(limits = c(.5, nn + .5), expand = c(0, 0)) +
  scale_fill_manual(values = c("0" = "#D6D6D6", "1" = unname(module_pal[1]), "2" = unname(module_pal[2]), "3" = unname(module_pal[3])), guide = "none") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(0, 4, 4, 4))
network <- (pt / bb) + plot_layout(heights = c(1, .055)) + plot_annotation(title = "A  Co-expression module discovery", subtitle = "50 patients | 5,000 variable genes", theme = theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.5)))
d <- subset(db$modules$current_associations, contrast == "B_C")
wr(d, "D_all_module_effects")
d$module <- factor(d$module, c("M03", "M02", "M01"))
d$label <- sprintf("adj. P %.3g", d$adj.P.Val)
effects <- ggplot(d, aes(logFC, module, colour = module)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#AAB4BE") +
  geom_segment(aes(x = CI.L, xend = CI.R, yend = module), linewidth = .6) +
  geom_point(size = 2) +
  geom_text(aes(x = CI.R, label = label), nudge_y = .24, hjust = 1, size = 2.3, colour = "#344453") +
  scale_color_manual(values = module_pal, guide = "none") +
  labs(title = "B  Associations with nonintubated ARDS", subtitle = "Nonintubated ARDS vs pneumonia without ARDS; all modules", x = "Module score difference (95% CI)", y = NULL)
sm <- w79$metadata$Sample[w79$discovery]
sm <- sm[order(w79$metadata$group[match(sm, w79$metadata$Sample)] == "B", w79$scores[sm, "M02"])]
genes <- w79$module_sets$M02
m <- t(scale(t(w79$expression[genes, sm])))
m <- m[hclust(dist(m))$order, , drop = FALSE]
genes <- rownames(m)
hh <- expand.grid(gene = genes, Sample = sm)
hh$z <- as.vector(m)
wr(hh, "E_expression")
hh$Sample <- factor(hh$Sample, sm)
hh$gene <- factor(hh$gene, genes)
# Recompute heatmap clustering from current 1,159-gene, 50-patient matrix.
stopifnot(length(sm) == 50, length(genes) == 1159, nrow(hh) == 50 * 1159)
ph <- ggplot(hh, aes(Sample, gene, fill = z)) +
  geom_raster() +
  scale_fill_gradient2(low = "#527D9A", mid = "white", high = "#B36562", limits = c(-2.5, 2.5), oob = scales::squish, name = "Gene z") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_text(size = 6))
top <- data.frame(Sample = factor(sm, sm), group = db$cohort$group[match(sm, db$cohort$Sample)])
wr(top, "E_sample_order")
bar <- ggplot(top, aes(Sample, 1, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = pal, labels = c(B = "Nonintubated ARDS (17)", C = "Pneumonia without ARDS (33)"), name = NULL) +
  theme_void() +
  theme(legend.position = "top", legend.text = element_text(size = 6))
heatmap <- (bar / ph) + plot_layout(heights = c(.10, 1)) + plot_annotation(title = "C  M02 expression structure", subtitle = "1,159 genes | discovery patients only", theme = theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.5)))
expanded <- readRDS("results/analysis/expanded_pathways/expanded_gsea_patient79.rds")
pick <- c("PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION", "PLATELET_AGGREGATION", "PLATELET_ACTIVATION", "RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2", "PLATELET_MORPHOGENESIS", "RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION", "MEGAKARYOCYTE_DEVELOPMENT", "MEGAKARYOCYTE_DIFFERENTIATION", "HEMOSTASIS", "CYTOPLASMIC_TRANSLATION", "EUKARYOTIC_TRANSLATION_INITIATION", "EUKARYOTIC_TRANSLATION_ELONGATION", "RIBOSOME", "COMPLEMENT_CASCADE", "COMPLEMENT_ACTIVATION", "COMPLEMENT_AND_COAGULATION_CASCADES", "NEUTROPHIL_DEGRANULATION", "LEUKOCYTE_MIGRATION")
x <- subset(expanded, contrast == "B_C" & pathway_label %in% pick)
x <- x[!x$pathway_label %in% c("NEUTROPHIL_DEGRANULATION", "LEUKOCYTE_MIGRATION") | x$pathway %in% c("REACTOME_NEUTROPHIL_DEGRANULATION", "GOBP_LEUKOCYTE_MIGRATION"), ]
x$label <- factor(wrap(pretty(x$pathway_label), 34), wrap(pretty(x$pathway_label[order(x$NES)]), 34))
levels(x$label) <- gsub("Runx1 Regulates Genes Involved in\\nMegakaryocyte Differentiation and\\nPlatelet Function", "RUNX1 regulation of\\nmegakaryocyte differentiation", levels(x$label))
x$category <- ifelse(grepl("NEUTROPHIL|LEUKOCYTE", x$pathway_label), "Immune / leukocyte", ifelse(grepl("COMPLEMENT", x$pathway_label), "Complement", ifelse(grepl("TRANSLATION|RIBOSOME", x$pathway_label), "Translation / ribosome", "Platelet / megakaryocyte / hemostasis")))
x$significant <- x$q_within_database < .05

h <- subset(db$gsea$primary, contrast %in% c("B_C", "A_C", "A_B"))
paths <- unique(h$pathway[h$padj < .05])
ref <- subset(h, contrast == "B_C" & pathway %in% paths)
lev <- ref$pathway[order(ref$NES)]
h <- subset(h, pathway %in% paths)
h$pathway <- factor(h$pathway, lev)
h$contrast <- factor(h$contrast, c("B_C", "A_C", "A_B"))
h$display_NES <- ifelse(h$padj < .05, h$NES, NA_real_)
wr(h, "A_Hallmark_contrasts")
pA <- ggplot(h, aes(contrast, pathway, fill = display_NES)) +
  geom_tile(colour = "white", linewidth = .4) +
  scale_fill_gradient2(low = "#527D9A", mid = "#FAFAFA", high = "#B36562", midpoint = 0, limits = c(-4, 4), na.value = "#E9ECEF", name = "NES") +
  scale_x_discrete(labels = cn, expand = c(0, 0)) +
  scale_y_discrete(labels = pretty, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6), legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = grid::unit(25, "mm"), barheight = grid::unit(2, "mm")))

b <- subset(x, category != "Complement")
b$lab <- wrap(pretty(b$pathway_label), 32)
b$lab[b$pathway_label == "RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION"] <- "RUNX1 regulation of\nmegakaryocyte differentiation"
b$lab[b$pathway_label == "PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION"] <- "Platelet activation, signaling\nand aggregation"
b$neglog10_q <- -log10(pmax(b$q_within_database, 1e-300))
b$program <- factor(b$category, levels = c("Platelet / megakaryocyte / hemostasis", "Immune / leukocyte", "Translation / ribosome"), labels = c("Platelet / hemostatic", "Immune", "Translation"))
# Actual facets with zero spacing; adjoining panel borders form thin separators.
b <- b[order(b$program, -b$NES, b$pathway), ]
b$lab <- factor(as.character(b$lab), levels = rev(as.character(b$lab)))
pBar <- ggplot(b, aes(neglog10_q, lab, fill = NES)) +
  geom_col(width = .55) +
  scale_fill_gradient2(low = "#527D9A", mid = "#FAFAFA", high = "#B36562", midpoint = 0, limits = c(-4, 4), name = "NES") +
  facet_grid(program ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, .06))) +
  scale_y_discrete(expand = expansion(add = .5)) +
  labs(x = expression(-log[10] ~ adjusted ~ italic(P)), y = NULL) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(axis.text.y = element_text(size = 5.5), axis.ticks.y = element_blank(), strip.text.y = element_text(size = 5.8, angle = -90, colour = "#465760", margin = margin(3, 3, 3, 3)), panel.spacing.y = grid::unit(0, "mm"), panel.border = element_rect(fill = NA, colour = "#222222", linewidth = .28), axis.line = element_blank(), strip.background = element_blank(), legend.position = "none")

pD <- effects + theme_classic(base_family = "Arial", base_size = 7) + theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_line(linewidth = .3, colour = "#64717A"))
wr(b, "B_displayed_expanded_pathways")
ps <- list(pA, pBar, network, pD, heatmap)
fig <- assemble(ps, design = "A#B\n###\nC#D\n###\nEEE", widths = c(1.08, .10, 1), heights = c(1.5, .025, .72, .025, .9))
finish(ps, fig, 3, 235, "V27 Figure 3. Full numerical results and post hoc display selection are exported as source data.")
