# 在项目根目录独立运行；输入为本次analysis目录中对应分析结果；输出为本图及源数据。
source("scripts/figures/common.R")
setup(4)
ann <- subset(db$modules$annotations, module == "M02")
wr(ann, "all_M02_annotations")
# Deterministic top-8 significance ranking with overlap-based redundancy removal.
select_terms <- function(x, n = 8) {
  x <- x[is.finite(x$q) & x$q < .05, ]
  x <- x[order(x$q, x$p, x$pathway), ]
  chosen <- integer()
  for (i in seq_len(nrow(x))) {
    g <- strsplit(x$genes[i], ";", fixed = TRUE)[[1]]
    redundant <- any(vapply(chosen, function(j) {
      h <- strsplit(x$genes[j], ";", fixed = TRUE)[[1]]
      length(intersect(g, h)) / length(union(g, h)) > .5
    }, logical(1)))
    if (!redundant) chosen <- c(chosen, i)
    if (length(chosen) >= n) break
  }
  x[chosen, , drop = FALSE]
}
add_metrics <- function(x) {
  module_size <- length(db$modules$sets$M02)
  background_size <- length(db$modules$background_genes)
  x$set_background <- vapply(x$pathway, function(path) length(intersect(db$gene_sets$sets[[path]], db$modules$background_genes)), integer(1))
  x$gene_ratio <- x$hits / module_size
  x$background_ratio <- x$set_background / background_size
  x$fold_enrichment <- x$gene_ratio / x$background_ratio
  x$coverage <- x$hits / x$set_background
  x$minus_log10_q <- -log10(pmax(x$q, 1e-300))
  x
}
# Functional categories are interpretive selections, not independent validation.
functional_all <- subset(ann, kind == "functional")
immune_refs <- c("REACTOME_INNATE_IMMUNE_SYSTEM", "REACTOME_NEUTROPHIL_DEGRANULATION", "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM", "GOBP_REGULATION_OF_LEUKOCYTE_MIGRATION")
platelet_terms <- select_terms(subset(functional_all, grepl("PLATELET|HEMOSTASIS|MEGAKARYOCYTE", pathway)), 4)
immune_terms <- functional_all[match(immune_refs, functional_all$pathway), ]
stopifnot(!anyNA(immune_terms$pathway), all(immune_terms$q < .05))
fun <- add_metrics(rbind(platelet_terms, immune_terms))
wr(functional_all, "A_complete_functional_context")
cell <- subset(ann, kind == "cell")
cell <- add_metrics(head(cell[order(cell$q, cell$p, cell$pathway), ], 6))
fun$label <- wrap(pretty(fun$pathway), 29)
cell$label <- wrap(pretty(cell$pathway), 29)
fun$label <- factor(fun$label, fun$label[order(fun$minus_log10_q)])
cell$label <- factor(cell$label, cell$label[order(cell$fold_enrichment)])
cell$stat <- sprintf("%d/%d genes", cell$hits, cell$set_background)
wr(fun, "A_selected_functions")
wr(cell, "B_selected_cell_references")
pA <- ggplot(fun, aes(minus_log10_q, label)) +
  geom_point(aes(size = hits, fill = fold_enrichment), shape = 21, colour = "white", stroke = .35, alpha = .95) +
  scale_size_continuous(range = c(3.2, 7.8), name = "Overlapping genes") +
  scale_fill_gradient(low = "#E2D9E0", high = unname(module_pal["M02"]), name = "Fold enrichment") +
  scale_x_continuous(expand = expansion(mult = c(.10, .08))) +
  labs(x = expression(-log[10] ~ adjusted ~ italic(P)), y = NULL) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "#30363B", linewidth = .35), axis.ticks = element_line(colour = "#30363B", linewidth = .35), panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 5.8), legend.position = "bottom", legend.box = "vertical") +
  guides(size = guide_legend(order = 1, title.position = "left", override.aes = list(fill = "#B05A92", colour = "white", alpha = .9)), fill = guide_colourbar(order = 2, title.position = "left", barwidth = grid::unit(22, "mm")))
pB <- ggplot(cell, aes(fold_enrichment, label, fill = coverage)) +
  geom_vline(xintercept = 1, colour = "#AAB4BE", linewidth = .3, linetype = 2) +
  geom_col(width = .58, colour = NA) +
  geom_text(aes(label = stat), hjust = 0, nudge_x = .06, size = 2.25, colour = "#344453") +
  scale_fill_gradient(low = "#DCE4E1", high = unname(module_pal["M02"]), limits = c(.4, 1), oob = scales::squish, breaks = c(.4, .6, .8, 1), labels = c("0.4", "0.6", "0.8", "1.0"), name = "Gene-set coverage") +
  scale_x_continuous(expand = expansion(mult = c(0, .34))) +
  labs(x = "Fold enrichment", y = NULL) +
  coord_cartesian(clip = "off") +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "#30363B", linewidth = .35), axis.ticks = element_line(colour = "#30363B", linewidth = .35), panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 5.8), legend.position = "bottom")
refs <- c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION", "GOBP_HEMOSTASIS", "GOBP_MEGAKARYOCYTE_DEVELOPMENT")
nm <- setNames(c("Platelet activation / aggregation", "Hemostasis", "Megakaryocyte development"), refs)
h <- subset(readRDS("results/analysis/expanded_pathways/expanded_gsea_patient79.rds"), contrast == "B_C" & pathway %in% refs)
le <- do.call(rbind, lapply(refs, function(path) {
  gg <- strsplit(h$leadingEdge[h$pathway == path], ";", fixed = TRUE)[[1]]
  data.frame(pathway = path, gene = gg, in_M02 = gg %in% db$modules$sets$M02)
}))
wr(le, "C_all_leading_edge_memberships")
od <- do.call(rbind, lapply(refs, function(path) {
  x <- subset(le, pathway == path)
  data.frame(pathway = path, category = c("M02", "Other genes"), count = c(sum(x$in_M02), sum(!x$in_M02)), total = nrow(x))
}))
od$pathway <- factor(od$pathway, rev(refs))
od$category <- factor(od$category, c("Other genes", "M02"))
wr(od, "C_overlap_counts")
lab <- subset(od, category == "M02")
lab$label <- paste0(lab$count, " / ", lab$total)
pC <- ggplot(od, aes(count, pathway, fill = category)) +
  geom_col(width = .5) +
  geom_text(data = lab, aes(x = total + 1, label = label), hjust = 0, size = 2.3) +
  scale_fill_manual(values = c("M02" = unname(module_pal["M02"]), "Other genes" = "#DCE2E5"), name = NULL, breaks = c("M02", "Other genes")) +
  scale_y_discrete(labels = nm) +
  scale_x_continuous(expand = expansion(mult = c(0, .20))) +
  labs(title = "C  Linking M02 to pathway enrichment", subtitle = "M02 / all leading-edge genes", x = "Leading-edge gene count", y = NULL)
sc <- db$modules$current_scores
sc$group <- factor(sc$group, c("D", "C", "B", "A"))
stopifnot(nrow(sc) == 113)
wr(sc, "D_all_group_scores")
pD <- ggplot(sc, aes(group, M02, colour = group)) +
  geom_boxplot(width = .5, outlier.shape = NA, linewidth = .4) +
  geom_point(position = position_jitter(width = .12, seed = 42), size = .9, alpha = .7) +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = c(D = "Healthy\nn=34", C = "Pneumonia\nwithout ARDS\nn=33", B = "Nonintubated ARDS\nn=17", A = "MV ARDS\nn=29")) +
  labs(title = "D  Extension across clinical states", x = NULL, y = "Projected M02 score") +
  theme(axis.text.x = element_text(size = 6))
d <- subset(db$modules$current_associations, module == "M02")
wr(d, "E_direct_contrasts")
d$label <- factor(cn[d$contrast], rev(cn))
d$stat <- sprintf("adj. P %.3g", d$adj.P.Val)
pE <- ggplot(d, aes(logFC, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#AAB4BE") +
  geom_segment(aes(x = CI.L, xend = CI.R, yend = label), colour = module_pal["M02"], linewidth = .6) +
  geom_point(colour = module_pal["M02"], size = 2) +
  geom_text(aes(x = CI.R, label = stat), nudge_y = .24, hjust = 1, size = 2.3) +
  labs(title = "E  Direct clinical-state comparisons", x = "M02 score difference (95% CI)", y = NULL)
fig <- assemble(list(pA, pB, pC, pD, pE), design = "AB\nCC\nDE", heights = c(1.35, .7, 1))
legend <- c(
  "Figure 4. Candidate platelet-related annotation and clinical-state variation of M02.",
  "A, Functional overrepresentation using the original 5,000 network genes as background. Horizontal position denotes −log10 BH-adjusted P, point size denotes overlapping gene count and colour denotes fold enrichment. Terms with 5–500 background-overlapping genes were tested. Four platelet/hemostatic terms selected by significance with hit-gene Jaccard redundancy filtering are displayed together with four explicit immune themes (innate immunity, neutrophil degranulation, cytokine signaling and leukocyte migration regulation). These interpretive categories were selected after analysis; all functional results are supplied. Enrichment measures module composition, not activation or direct cell interaction. Complete results are supplied.",
  "B, The six reference cell-type sets with smallest adjusted P among the 54 healthy lung single-cell reference sets, subject to the original size filter. Bar length denotes fold enrichment and colour denotes the fraction of each reference set represented in M02. Labels show overlapping genes/reference-set genes. Cell-reference overlap supports annotation rather than exclusive tissue or cell provenance.",
  "C, Descriptive overlap of the patient-only platelet activation/aggregation, hemostasis and megakaryocyte-development leading edges with original M02. Labels report M02 members/all leading-edge members; full memberships are supplied.",
  "D, Fixed discovery-loading M02 scores across 115 profiles. Points denote profiles; boxes depict median/IQR and whiskers 1.5 IQR. The primary display retains the pre-intubation profile of the cross-state pair. Healthy controls are a reference, not a temporal baseline.",
  "E, Direct group contrasts with 95% confidence intervals and BH adjustment across three modules per contrast. Clinical-state differences do not establish a ventilation effect. All source data are supplied."
)
finish(list(pA, pB, pC, pD, pE), fig, 4, 238, legend, list(c(150, 125), c(145, 110), c(183, 65), c(130, 95), c(130, 95)))
