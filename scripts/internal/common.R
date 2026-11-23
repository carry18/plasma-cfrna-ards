source("scripts/config.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})
root <- input_dir
refs <- readRDS(file.path(input_dir, "reference_sets.rds"))
sets <- refs$hallmark
z <- function(x) as.numeric(scale(x))
score_expr <- function(e, gs) {
  sdg <- apply(e, 1, sd)
  e <- e[is.finite(sdg) & sdg > 0, , drop = FALSE]
  gz <- t(scale(t(e)))
  scores <- sapply(gs, function(g) colMeans(gz[intersect(g, rownames(gz)), , drop = FALSE]))
  list(scores = scores, z = gz, coverage = sapply(gs, function(g) length(intersect(g, rownames(gz)))))
}
