#!/usr/bin/env Rscript

source("scripts/internal/common.R")
suppressPackageStartupMessages(library(limma))

out <- "results/analysis/wgcna"
checkpoint <- readRDS(file.path(out, "analysis_checkpoint_before_external.rds"))
external <- readRDS(file.path(root, "public_context_processed.rds"))$GSE66890
focal_module <- checkpoint$focal_module
focal_genes <- checkpoint$module_sets[[focal_module]]

external_matrix <- external$expr[intersect(focal_genes, rownames(external$expr)), , drop = FALSE]
external_complete <- rowSums(is.finite(external_matrix)) == ncol(external_matrix)
external_variable <- apply(external_matrix, 1, sd, na.rm = TRUE) > 0
external_genes <- rownames(external_matrix)[external_complete & external_variable]
external_score <- colMeans(t(scale(t(external_matrix[external_genes, , drop = FALSE]))))
external_score <- as.numeric(scale(external_score))
external_design <- model.matrix(~ ARDS + age + sex + batch, external$metadata)
external_fit <- lmFit(matrix(external_score, nrow = 1), external_design)
external_index <- match("ARDS", colnames(external_design))
external_estimate <- external_fit$coefficients[1, external_index]
external_se <- external_fit$stdev.unscaled[1, external_index] * external_fit$sigma[1]
external_df <- external_fit$df.residual[1]
external_t <- external_estimate / external_se
external_p <- 2 * pt(abs(external_t), df = external_df, lower.tail = FALSE)
external_critical <- qt(0.975, df = external_df)
external_effect <- data.frame(
  module = focal_module,
  measured_genes = length(external_genes),
  estimate = external_estimate,
  lower = external_estimate - external_critical * external_se,
  upper = external_estimate + external_critical * external_se,
  P = external_p
)
write.csv(external_effect, file.path(out, "GSE66890_focal_score_effect.csv"), row.names = FALSE)

patient_md <- droplevels(checkpoint$metadata[
  checkpoint$metadata$group %in% c("A", "B", "C"),
])
patient_expression <- checkpoint$expression[, rownames(checkpoint$scores) %in% patient_md$Sample,
  drop = FALSE
]
patient_expression <- patient_expression[, patient_md$Sample, drop = FALSE]
shared <- intersect(focal_genes, intersect(rownames(patient_expression), rownames(external$expr)))
shared <- shared[rowSums(is.finite(external$expr[shared, , drop = FALSE])) == ncol(external$expr)]
internal_residual <- t(qr.resid(
  qr(model.matrix(~group, patient_md)),
  t(patient_expression[shared, , drop = FALSE])
))
external_residual <- t(qr.resid(
  qr(external_design),
  t(external$expr[shared, , drop = FALSE])
))
shared <- shared[apply(internal_residual[shared, , drop = FALSE], 1, sd) > 0 &
  apply(external_residual[shared, , drop = FALSE], 1, sd) > 0]
internal_z <- t(scale(t(internal_residual[shared, , drop = FALSE])))
external_z <- t(scale(t(external_residual[shared, , drop = FALSE])))
kme <- data.frame(
  gene = shared,
  internal_kME = as.numeric(cor(t(internal_z), colMeans(internal_z))),
  external_kME = as.numeric(cor(t(external_z), colMeans(external_z)))
)
write.csv(kme, file.path(out, "GSE66890_gene_module_correspondence.csv"), row.names = FALSE)

checkpoint$external <- list(effect = external_effect, kME = kme)
saveRDS(checkpoint, file.path(out, "analysis.rds"))
cat("Focal module:", focal_module, "\n")
print(external_effect)
cat("External kME Pearson:", cor(kme$internal_kME, kme$external_kME), "\n")
cat("External kME Spearman:", cor(kme$internal_kME, kme$external_kME,
  method = "spearman"
), "\n")
