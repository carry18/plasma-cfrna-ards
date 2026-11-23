# 在项目根目录独立运行；输入为本次analysis目录中对应分析结果；输出为本图及源数据。
source("scripts/figures/common.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})
setup(2)
g <- subset(db$differential_expression, dataset == "internal_cfRNA" & contrast == "ABC_D" & model == "base")
stopifnot(!anyDuplicated(g$gene))
lfc_cutoff <- 1
g$status <- ifelse(g$padj < .05 & g$logFC >= lfc_cutoff, "Higher in patients", ifelse(g$padj < .05 & g$logFC <= -lfc_cutoff, "Lower in patients", "Other genes"))
g$status <- factor(g$status, c("Higher in patients", "Lower in patients", "Other genes"))
g$lp <- -log10(pmax(g$padj, 1e-300))
wr(g, "A_all_patient_vs_healthy_genes")
wr(data.frame(padj_cutoff = .05, absolute_log2FC_cutoff = lfc_cutoff, higher = sum(g$status == "Higher in patients"), lower = sum(g$status == "Lower in patients"), other = sum(g$status == "Other genes")), "A_thresholds")
pA <- ggplot(g, aes(logFC, lp)) +
  geom_point(aes(colour = status), size = .5, alpha = .65) +
  geom_hline(yintercept = -log10(.05), linetype = 2, colour = "#8B999F", linewidth = .3) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, colour = "#8B999F", linewidth = .3) +
  scale_colour_manual(values = c("Higher in patients" = "#B36562", "Lower in patients" = "#527D9A", "Other genes" = "#CCD2D6"), name = NULL, drop = FALSE) +
  labs(x = "log2 fold change", y = "−log10 adjusted P") +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "#829099", linewidth = .35), legend.position = "bottom", legend.box = "horizontal") +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1.5, alpha = 1)))
# Top 60 absolute moderated-t genes among adjusted-significant ABC_D genes.
sig <- g[g$padj < .05, ]
sig <- sig[order(-abs(sig$statistic), sig$gene), ]
genes <- head(sig$gene, 60)
wr(sig[sig$gene %in% genes, ], "B_selected60_genes")
md <- db$cohort
md$group <- factor(md$group, c("A", "B", "C", "D"))
md$Batch <- factor(md$Batch)
e <- if (!is.null(w79)) w79$expression else removeBatchEffect(
  a79$descriptive$expression, batch = md$Batch, design = model.matrix(~group, md)
)

sm <- md$Sample[order(match(md$group, c("D", "C", "B", "A")), md$Sample)]
m <- t(scale(t(e[genes, sm])))
m <- m[hclust(dist(m))$order, , drop = FALSE]
x <- expand.grid(gene = rownames(m), Sample = sm)
x$z <- as.vector(m)
x$gene <- factor(x$gene, rev(rownames(m)))
x$Sample <- factor(x$Sample, sm)
wr(x, "B_expression60")
top <- data.frame(Sample = factor(sm, sm), group = md$group[match(sm, md$Sample)])
wr(top, "B_sample_order")
heat <- ggplot(x, aes(Sample, gene, fill = z)) +
  geom_raster() +
  scale_fill_gradient2(low = "#527D9A", mid = "white", high = "#B36562", limits = c(-2.5, 2.5), oob = scales::squish, name = "Gene z") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_text(size = 6))
bar <- ggplot(top, aes(Sample, 1, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = pal, labels = c(A = "MV ARDS (29)", B = "Nonintubated ARDS (17)", C = "Pneumonia without ARDS (33)", D = "Healthy (34)"), name = NULL) +
  theme_void(base_family = "Arial") +
  theme(legend.position = "top", legend.text = element_text(size = 5.7)) +
  guides(fill = guide_legend(nrow = 2))
pB <- (bar / heat) + plot_layout(heights = c(.10, 1))
de3 <- subset(db$differential_expression, dataset == "internal_cfRNA" & model == "base" & contrast %in% c("A_D", "B_D", "C_D"))
wide <- reshape(de3[, c("gene", "contrast", "logFC", "padj")], idvar = "gene", timevar = "contrast", direction = "wide")
wide$A <- wide$padj.A_D < .05
wide$B <- wide$padj.B_D < .05
wide$C <- wide$padj.C_D < .05
wide$key <- paste0(ifelse(wide$A, "A", ""), ifelse(wide$B, "B", ""), ifelse(wide$C, "C", ""))
wide$direction <- ifelse(wide$logFC.A_D > 0 & wide$logFC.B_D > 0 & wide$logFC.C_D > 0, "Higher in patients", ifelse(wide$logFC.A_D < 0 & wide$logFC.B_D < 0 & wide$logFC.C_D < 0, "Lower in patients", "Discordant"))
wide_sig <- wide[wide$key != "", ]
ic <- as.data.frame(table(wide_sig$key), stringsAsFactors = FALSE)
names(ic) <- c("key", "n")
ic <- ic[order(-ic$n), ]
ic$key <- factor(ic$key, levels = as.character(ic$key))
wr(ic, "C_exact_intersection_counts")
shared <- wide_sig[wide_sig$key == "ABC", ]
sd <- as.data.frame(table(factor(shared$direction, levels = c("Higher in patients", "Lower in patients", "Discordant"))), stringsAsFactors = FALSE)
names(sd) <- c("direction", "n")
sd <- sd[sd$n > 0, ]
wr(sd, "C_shared_gene_directions")
cors <- data.frame(comparison = c("MV ARDS vs nonintubated ARDS", "MV ARDS vs pneumonia without ARDS", "Nonintubated ARDS vs pneumonia without ARDS"), pearson_r = c(cor(wide$logFC.A_D, wide$logFC.B_D), cor(wide$logFC.A_D, wide$logFC.C_D), cor(wide$logFC.B_D, wide$logFC.C_D)))
wr(cors, "C_all_gene_logFC_correlations")
wr(wide, "C_gene_level_membership_and_effects")
venn_sets <- list("MV ARDS" = wide$gene[wide$A], "Nonintubated ARDS" = wide$gene[wide$B], "Pneumonia without ARDS" = wide$gene[wide$C])
venn_glist <- suppressMessages(VennDiagram::venn.diagram(x = venn_sets, filename = NULL, disable.logging = TRUE, euler.d = FALSE, scaled = FALSE, fill = unname(pal[c("A", "B", "C")]), col = unname(pal[c("A", "B", "C")]), alpha = c(.34, .34, .34), lwd = .82, cex = .68, fontfamily = "Arial", fontface = "plain", cat.col = unname(pal[c("A", "B", "C")]), cat.cex = .62, cat.fontfamily = "Arial", cat.fontface = "plain", cat.pos = c(-22, 22, 180), cat.dist = c(.055, .055, .045), margin = .07))
pC_venn <- wrap_elements(full = grid::gTree(children = venn_glist, vp = grid::viewport(width = grid::unit(.98, "snpc"), height = grid::unit(.98, "snpc")))) + theme(plot.background = element_rect(fill = "white", colour = NA), panel.border = element_blank(), plot.margin = margin(3, 3, 3, 3))
sd_plot <- data.frame(direction = factor(c("Consistently lower vs healthy", "Consistently higher vs healthy"), levels = c("Consistently lower vs healthy", "Consistently higher vs healthy")), n = c(sum(shared$direction == "Lower in patients"), sum(shared$direction == "Higher in patients")), label = c("Lower", "Higher"))
wr(sd_plot, "C_shared_direction_summary")
pC <- pC_venn / plot_spacer() + plot_layout(heights = c(1, .001))
# Direction-aware ORA of the shared genes shared by all three disease-versus-healthy contrasts.
pathway_ora <- function(query_genes, direction, collection, prefix) {
  bg <- unique(wide$gene)
  sets <- db$gene_sets$sets
  sn <- names(sets)[startsWith(names(sets), prefix)]
  z <- do.call(rbind, lapply(sn, function(pathway) {
    set_genes <- intersect(unique(sets[[pathway]]), bg)
    n_set <- length(set_genes)
    if (n_set < 10L || n_set > 500L) {
      return(NULL)
    }
    hit <- intersect(query_genes, set_genes)
    expected <- length(query_genes) * n_set / length(bg)
    data.frame(collection = collection, direction = direction, pathway = pathway, query_size = length(query_genes), background_size = length(bg), set_size = n_set, overlap = length(hit), expected = expected, fold_enrichment = length(hit) / expected, p = phyper(length(hit) - 1L, n_set, length(bg) - n_set, length(query_genes), lower.tail = FALSE), genes = paste(sort(hit), collapse = ";"), stringsAsFactors = FALSE)
  }))
  z$padj <- p.adjust(z$p, "BH")
  z[order(z$padj, z$p), ]
}
reactome_all <- rbind(
  pathway_ora(shared$gene[shared$direction == "Higher in patients"], "Consistently higher", "Reactome", "REACTOME_"),
  pathway_ora(shared$gene[shared$direction == "Lower in patients"], "Consistently lower", "Reactome", "REACTOME_")
)
hallmark_all <- rbind(
  pathway_ora(shared$gene[shared$direction == "Higher in patients"], "Consistently higher", "Hallmark", "HALLMARK_"),
  pathway_ora(shared$gene[shared$direction == "Lower in patients"], "Consistently lower", "Hallmark", "HALLMARK_")
)
wr(reactome_all, "D_all_reactome_ORA")
wr(hallmark_all, "D_all_hallmark_ORA")
d_keep <- c("HALLMARK_ANGIOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "REACTOME_FIBRONECTIN_MATRIX_FORMATION", "REACTOME_MITOCHONDRIAL_TRANSLATION", "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT")
d_direction <- c(rep("Consistently higher", 3), rep("Consistently lower", 2))
d_pool <- rbind(hallmark_all, reactome_all)
d_show <- do.call(rbind, lapply(seq_along(d_keep), function(i) d_pool[d_pool$pathway == d_keep[i] & d_pool$direction == d_direction[i], , drop = FALSE]))
stopifnot(nrow(d_show) == length(d_keep), all(is.finite(d_show$padj)))
d_names <- c("Angiogenesis", "Epithelial–mesenchymal transition", "Fibronectin matrix formation", "Mitochondrial translation", "Aerobic respiration and\nrespiratory electron transport")
d_show$label <- factor(d_names, levels = rev(d_names))
d_show$neglog10_padj <- -log10(pmax(d_show$padj, 1e-300))
wr(d_show, "D_displayed_curated_pathways")
pD <- ggplot(d_show, aes(neglog10_padj, label, size = fold_enrichment, colour = direction)) +
  geom_vline(xintercept = -log10(.05), linetype = 2, colour = "#8B999F", linewidth = .3) +
  geom_point(alpha = .65) +
  scale_x_continuous(expand = expansion(mult = c(.04, .16))) +
  scale_colour_manual(values = c("Consistently higher" = "#B36562", "Consistently lower" = "#527D9A"), breaks = c("Consistently higher", "Consistently lower"), labels = c("Higher in patients", "Lower in patients"), name = NULL) +
  scale_size_continuous(range = c(3, 7), breaks = c(3, 7, 11), name = "Fold enrichment") +
  labs(x = "−log10 adjusted P", y = NULL) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "#829099", linewidth = .35), axis.text.y = element_text(size = 5.55), panel.grid.major.y = element_line(colour = "#EDF0F2", linewidth = .25), legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 5.3), legend.title = element_text(size = 5.5), legend.spacing.y = grid::unit(.5, "mm"), legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(order = 1, nrow = 1, byrow = TRUE, override.aes = list(size = 3, alpha = 1)), size = guide_legend(order = 2, nrow = 1))
ab <- subset(db$gsea$primary, contrast == "AB_C" & model == "base")
wr(ab, "E_all50_aggregate_ARDS")
aa <- subset(ab, padj < .05)
wr(aa, "E_all_significant_aggregate")
dd <- aa
dd$lp <- -log10(pmax(dd$padj, 1e-300))
dd$pathway <- factor(dd$pathway, dd$pathway[order(dd$NES)])
d_lab <- setNames(pretty(levels(dd$pathway)), levels(dd$pathway))
# All significant Hallmarks shown without privileging three pathways.
wr(dd, "E_displayed_aggregate")
pE <- ggplot(dd, aes(NES, pathway)) +
  geom_vline(xintercept = 0, colour = "#829099", linewidth = .35) +
  geom_col(aes(fill = NES), width = .68) +
  scale_fill_gradient2(low = "#527D9A", mid = "#ECEFF1", high = "#B36562", midpoint = 0, name = "NES") +
  scale_y_discrete(labels = d_lab) +
  labs(x = "Normalized enrichment score", y = NULL) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "#829099", linewidth = .35), axis.text.y = element_text(size = 5.15), panel.grid.major.y = element_blank(), legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = grid::unit(24, "mm")))
fig <- assemble(list(pA, pB, pC, pD, pE), design = "AB\nCE\nDE", widths = c(.92, 1.08), heights = c(.74, .76, 1.04))
legend <- c(
  "Figure 2. Disease-associated plasma cfRNA remodeling and aggregate ARDS pathway patterns.",
  "A, Equal-weighted patient-group mean versus healthy controls in a group- and batch-adjusted descriptive model of 113 profiles (79 patients and 34 healthy controls). All tested genes are shown. Lines mark BH-adjusted P=0.05 and absolute log2 fold change=1.",
  "B, The 60 adjusted-significant genes with largest absolute moderated t in A. Batch-adjusted expression preserves group effects; gene-z colour saturates at ±2.5. All 113 profiles are retained. This is a descriptive display of the same cohort.",
  sprintf("C, Fixed-circle Venn diagram of BH-significant genes in the three patient-versus-healthy contrasts. %d genes are shared: %d consistently higher, %d consistently lower and %d discordant. Circle area is not proportional to set size.", nrow(shared), sum(shared$direction == "Higher in patients"), sum(shared$direction == "Lower in patients"), sum(shared$direction == "Discordant")),
  "D, Direction-aware Hallmark and Reactome overrepresentation of shared genes using the tested-gene universe. BH adjustment is separate within collection and direction. Position shows −log10 adjusted P, size fold enrichment and colour direction. Five displayed themes are manually selected after analysis; complete results accompany the figure.",
  sprintf("E, All %d BH-significant Hallmark sets in the equal-weighted overall ARDS-versus-pneumonia without ARDS contrast, fitted to 79 patients with batch adjustment. Significance is corrected across 50 Hallmarks. Aggregate effects motivate subsequent state-specific analysis.", nrow(dd))
)
finish(list(pA, pB, pC, pD, pE), fig, 2, 245, legend, list(c(130, 100), c(150, 110), c(145, 115), c(150, 115), c(150, 190)))
