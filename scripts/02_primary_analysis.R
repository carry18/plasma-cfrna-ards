# 输入：01生成的counts、analysis_samples.csv及基因集。
# 输出：差异表达、标准化表达、通路富集、固定评分及分析样本信息。
# 在项目根目录运行；前置步骤：01。
#!/usr/bin/env Rscript

# Analysis inclusion is defined by the final upload sample manifest.

source("scripts/internal/common.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
})

set.seed(20261016)
out <- "results/analysis"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "expanded_pathways"), recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(root, "internal_processed.rds"))
inc <- read.csv("results/inputs/analysis_samples.csv",
  check.names = FALSE
)
inc$group <- factor(inc$group, levels = c("A", "B", "C", "D"))
inc$Batch <- factor(inc$Batch)

main_md <- droplevels(inc[inc$group %in% c("A", "B", "C"), ])
descr_md <- droplevels(inc)
stopifnot(nrow(main_md) == 79L, nrow(descr_md) == 113L)
stopifnot(identical(as.integer(table(main_md$group)), c(29L, 17L, 33L)))

main_comparisons <- c(
  AB_C = "(groupA+groupB)/2-groupC",
  B_C = "groupB-groupC",
  A_C = "groupA-groupC",
  A_B = "groupA-groupB"
)
descr_comparisons <- c(
  A_D = "groupA-groupD",
  B_D = "groupB-groupD",
  C_D = "groupC-groupD",
  ABC_D = "(groupA+groupB+groupC)/3-groupD"
)

run_model <- function(md, comparisons, cohort) {
  ct <- d$counts[, as.character(md$Sample), drop = FALSE]
  dg <- DGEList(ct)
  dg <- dg[filterByExpr(dg, group = md$group), , keep.lib.sizes = FALSE]
  dg <- calcNormFactors(dg)
  design <- model.matrix(~ 0 + group + Batch, md)
  stopifnot(qr(design)$rank == ncol(design))
  v <- voom(dg, design, plot = FALSE)
  contrast_matrix <- makeContrasts(contrasts = comparisons, levels = design)
  colnames(contrast_matrix) <- names(comparisons)
  fit <- eBayes(contrasts.fit(lmFit(v, design), contrast_matrix), robust = TRUE)

  genes <- list()
  enrichment <- list()
  summary <- list()
  for (contrast in names(comparisons)) {
    tab <- topTable(fit, coef = contrast, number = Inf, sort.by = "none")
    tab$gene <- rownames(tab)
    tab$contrast <- contrast
    tab$cohort <- cohort
    genes[[paste(cohort, contrast, sep = "|")]] <- tab
    summary[[contrast]] <- data.frame(
      cohort = cohort,
      contrast = contrast,
      tested_genes = nrow(tab),
      nominal_p_lt_0_05 = sum(tab$P.Value < 0.05),
      adjusted_p_lt_0_05 = sum(tab$adj.P.Val < 0.05),
      stringsAsFactors = FALSE
    )

    ranks <- sort(setNames(tab$t, tab$gene), decreasing = TRUE)
    set.seed(20261016)
    fg <- as.data.frame(fgseaMultilevel(
      sets, ranks,
      minSize = 10, maxSize = 500, eps = 1e-10,
      nPermSimple = 10000, BPPARAM = BiocParallel::SerialParam()
    ))
    fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = "")
    fg$contrast <- contrast
    fg$cohort <- cohort
    enrichment[[contrast]] <- fg
    message(cohort, " ", contrast, " complete")
  }

  list(
    metadata = md,
    genes = genes,
    enrichment = do.call(rbind, enrichment),
    summary = do.call(rbind, summary),
    expression = v$E,
    design = design
  )
}

main <- run_model(main_md, main_comparisons, "unique_patient79")
descriptive <- run_model(descr_md, descr_comparisons, "unique_profile113")
gene_results <- c(main$genes, descriptive$genes)
hallmark_results <- rbind(main$enrichment, descriptive$enrichment)

saveRDS(gene_results, file.path(out, "gene_results.rds"))
saveRDS(
  list(main = main, descriptive = descriptive, sets = sets),
  file.path(out, "analysis_inputs.rds")
)
write.csv(hallmark_results, file.path(out, "hallmark_all_contrasts.csv"), row.names = FALSE)
write.csv(rbind(main$summary, descriptive$summary),
  file.path(out, "gene_summary.csv"),
  row.names = FALSE
)
write.csv(main_md, file.path(out, "included_samples_patients79.csv"), row.names = FALSE)
write.csv(descr_md, file.path(out, "included_samples_descriptive113.csv"), row.names = FALSE)
write.csv(as.data.frame(table(main_md$group, main_md$Batch)),
  file.path(out, "group_batch.csv"),
  row.names = FALSE
)

# Expanded Reactome/GO/KEGG enrichment for patient contrasts.
expanded_sets <- readRDS("results/inputs/expanded_sets.rds")
definitions <- read.csv("results/inputs/pathway_definitions.csv")
universe <- unique(main$genes[[1]]$gene)
expanded_sets <- lapply(expanded_sets, intersect, universe)
expanded_sets <- expanded_sets[lengths(expanded_sets) >= 10 & lengths(expanded_sets) <= 500]
database <- setNames(definitions$database, definitions$pathway)[names(expanded_sets)]
stopifnot(all(!is.na(database)))

expanded_results <- list()
for (contrast in names(main_comparisons)) {
  tab <- main$genes[[paste("unique_patient79", contrast, sep = "|")]]
  ranks <- sort(setNames(tab$t, tab$gene), decreasing = TRUE)
  set.seed(20261016)
  fg <- as.data.frame(fgseaMultilevel(
    expanded_sets, ranks,
    minSize = 10, maxSize = 500, eps = 1e-4,
    nPermSimple = 1000, BPPARAM = BiocParallel::SerialParam()
  ))
  fg$database <- database[fg$pathway]
  fg$q_within_database <- ave(fg$pval, fg$database,
    FUN = function(p) p.adjust(p, method = "BH")
  )
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = "")
  fg$contrast <- contrast
  fg$pathway_label <- sub("^(GOBP|REACTOME|KEGG|KEGG_MEDICUS)_", "", fg$pathway)
  expanded_results[[contrast]] <- fg
  message("Expanded ", contrast, " complete")
}
expanded_results <- do.call(rbind, expanded_results)
write.csv(expanded_results,
  file.path(out, "expanded_pathways", "expanded_gsea_patient79.csv"),
  row.names = FALSE
)
saveRDS(
  expanded_results,
  file.path(out, "expanded_pathways", "expanded_gsea_patient79.rds")
)

# Fixed mean-z pathway scores in the unique-patient cohort.
hallmark_def <- read.csv("results/inputs/hallmark_definition.csv")
hallmark_sets <- lapply(split(hallmark_def$gene_symbol, hallmark_def$gs_name), unique)
pick <- function(collection, ids) setNames(lapply(ids, function(id) collection[[id]]), ids)
targets <- c(
  pick(hallmark_sets, c(
    "HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_COAGULATION"
  )),
  pick(expanded_sets, c(
    "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
    "GOBP_PLATELET_ACTIVATION", "GOBP_PLATELET_AGGREGATION",
    "REACTOME_HEMOSTASIS", "GOBP_HEMOSTASIS",
    "GOBP_MEGAKARYOCYTE_DEVELOPMENT", "GOBP_MEGAKARYOCYTE_DIFFERENTIATION",
    "GOBP_HEMOPOIESIS", "GOBP_ERYTHROCYTE_DEVELOPMENT",
    "GOBP_CYTOPLASMIC_TRANSLATION", "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"
  ))
)

zmat <- t(scale(t(main$expression)))
score_matrix <- sapply(targets, function(gene_set) {
  genes <- intersect(gene_set, rownames(zmat))
  colMeans(zmat[genes, , drop = FALSE])
})
score_matrix <- scale(score_matrix)
coverage <- sapply(targets, function(gene_set) length(intersect(gene_set, rownames(zmat))))
score_contrasts <- c(
  B_C = "groupB-groupC", A_C = "groupA-groupC", A_B = "groupA-groupB",
  AB_C = "(groupA+groupB)/2-groupC"
)
contrast_matrix <- makeContrasts(contrasts = score_contrasts, levels = main$design)
colnames(contrast_matrix) <- names(score_contrasts)
score_fit <- eBayes(contrasts.fit(lmFit(t(score_matrix), main$design), contrast_matrix),
  robust = TRUE
)
score_results <- do.call(rbind, lapply(names(score_contrasts), function(contrast) {
  tab <- topTable(score_fit, coef = contrast, number = Inf, sort.by = "none", confint = TRUE)
  data.frame(
    contrast = contrast,
    score = rownames(tab),
    estimate = tab$logFC,
    lower = tab$CI.L,
    upper = tab$CI.R,
    P = tab$P.Value,
    adjusted_P = tab$adj.P.Val,
    measured_genes = coverage[rownames(tab)],
    row.names = NULL
  )
}))
write.csv(score_results, file.path(out, "fixed_score_effects.csv"), row.names = FALSE)
write.csv(
  data.frame(
    Sample = rownames(score_matrix), group = main_md$group,
    Batch = main_md$Batch, score_matrix, check.names = FALSE
  ),
  file.path(out, "fixed_score_values.csv"),
  row.names = FALSE
)

hallmark_summary <- aggregate(
  hallmark_results$padj < 0.05,
  list(cohort = hallmark_results$cohort, contrast = hallmark_results$contrast),
  sum,
  na.rm = TRUE
)
names(hallmark_summary)[3] <- "significant_hallmark_sets"
write.csv(hallmark_summary, file.path(out, "hallmark_significant_counts.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("\nGroup counts\n")
print(table(main_md$group, main_md$Batch))
cat("\nSignificant Hallmark sets\n")
print(hallmark_summary[hallmark_summary$cohort == "unique_patient79", ])
cat("DONE\n")
