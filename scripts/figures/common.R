# Shared style/export helpers only; each FigureN.R can run independently from project root.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})
base_dir <- "results/figures"
db <- list(gsea = list(), modules = list(), gene_sets = list())
pal <- c(D = "#8B999F", C = "#527D9A", B = "#93749E", A = "#B38A56")
module_pal <- c(M01 = "#3F8F83", M02 = "#B05A92", M03 = "#8A9147")
cn <- c(B_C = "Nonintubated ARDS\nvs pneumonia\nwithout ARDS", A_C = "MV ARDS\nvs pneumonia\nwithout ARDS", A_B = "MV vs nonintubated\nARDS")
th <- theme_bw(base_family = "Arial", base_size = 7) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#EDF0F2", linewidth = .2), panel.border = element_rect(colour = "#829099", linewidth = .35), plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.3, colour = "#52616B"), legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_text(size = 6.5), plot.margin = margin(6, 6, 6, 6), strip.background = element_rect(fill = "#F4F6F7", colour = NA), strip.text = element_text(size = 7))
theme_set(th)
options(device = function(...) ragg::agg_png(tempfile(fileext = ".png")))
setup <- function(n) {
  out <<- file.path(base_dir, paste0("Figure", n))
  dir.create(file.path(out, "source_data"), recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = normalizePath(out))
}
wr <- function(x, n) write.csv(x, file.path(out, "source_data", paste0(n, ".csv")), row.names = FALSE)
pretty <- function(v) {
  v <- tools::toTitleCase(tolower(gsub("_", " ", sub("^(HALLMARK_|GOBP_|REACTOME_|TRAVAGLINI_LUNG_)", "", v))))
  fix <- c(Tnfa = "TNF-alpha", Nfkb = "NF-kB", Pi3k = "PI3K", Akt = "AKT", Mtorc1 = "mTORC1", Mtor = "mTOR", Myc = "MYC", G2m = "G2M", E2f = "E2F", Dna = "DNA", Il6 = "IL6", Jak = "JAK", Stat3 = "STAT3", Uv = "UV", Tgf = "TGF", Kras = "KRAS")
  for (k in names(fix)) v <- gsub(k, fix[[k]], v, fixed = TRUE)
  v
}
wrap <- function(v, width = 32) vapply(v, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
# Publication policy: descriptive titles, subtitles and captions belong in legends.
clean_panel <- function(p) {
  if (inherits(p, "patchwork")) {
    if (length(p$patches$plots)) p$patches$plots <- lapply(p$patches$plots, clean_panel)
    p$patches$annotation$title <- NULL
    p$patches$annotation$subtitle <- NULL
    p$patches$annotation$caption <- NULL
  }
  p + labs(title = NULL, subtitle = NULL, caption = NULL) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())
}
flat <- function(p) wrap_elements(full = if (inherits(p, "patchwork")) patchworkGrob(p) else ggplotGrob(p))
letter_panel <- function(p, label) {
  flat(clean_panel(p)) + labs(tag = label) + theme(plot.tag = element_text(family = "Arial", face = "bold", size = 9), plot.tag.position = c(0, 1), plot.margin = margin(5, 5, 5, 5))
}
assemble <- function(ps, design, widths = NULL, heights = NULL) wrap_plots(lapply(seq_along(ps), function(i) letter_panel(ps[[i]], LETTERS[i])), design = design, widths = widths, heights = heights)
export <- function(p, name, w = 183, h = 180, tiff = TRUE) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  svglite(file.path(out, paste0(name, ".svg")), width = w / 25.4, height = h / 25.4)
  print(p)
  dev.off()
  cairo_pdf(file.path(out, paste0(name, ".pdf")), width = w / 25.4, height = h / 25.4, family = "Arial")
  print(p)
  dev.off()
  if (tiff) {
    agg_tiff(file.path(out, paste0(name, ".tiff")), width = w, height = h, units = "mm", res = 600, compression = "lzw")
    print(p)
    dev.off()
  }
}
finish <- function(ps, fig, n, h, legend, panel_dims = NULL) {
  ps <- lapply(ps, clean_panel)
  export(fig, paste0("Figure", n), h = h)
  for (i in seq_along(ps)) {
    dim <- if (is.null(panel_dims)) c(130, 100) else panel_dims[[i]]
    export(letter_panel(ps[[i]], LETTERS[i]), paste0("Figure", n, LETTERS[i]), w = dim[1], h = dim[2], tiff = FALSE)
  }
  saveRDS(ps, file.path(out, "plots.rds"))
  writeLines(legend, file.path(out, "legend.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}


a79 <- readRDS("results/analysis/analysis_inputs.rds")
network_path <- "results/analysis/wgcna/analysis.rds"
w79 <- if (file.exists(network_path)) readRDS(network_path) else NULL
gnew <- do.call(rbind, readRDS("results/analysis/gene_results.rds"))
gnew$padj <- gnew$adj.P.Val
gnew$statistic <- gnew$t
gnew$dataset <- "internal_cfRNA"
gnew$model <- "base"
db$differential_expression <- gnew
h79 <- read.csv("results/analysis/hallmark_all_contrasts.csv")
db$gsea$primary <- subset(h79, cohort == "unique_patient79")
db$gsea$primary$model <- "base"
db$cohort <- a79$descriptive$metadata
if (!is.null(w79)) {
db$cohort <- w79$metadata
db$modules$sets <- w79$module_sets
db$modules$background_genes <- w79$network_genes
db$modules$current_scores <- cbind(w79$metadata, as.data.frame(w79$scores))
z <- w79$effects
z$logFC <- z$estimate
z$CI.L <- z$lower
z$CI.R <- z$upper
z$adj.P.Val <- z$adjusted_P
db$modules$current_associations <- subset(z, contrast != "AB_C")
z <- w79$annotations
z$q <- z$adjusted_P
z$p <- z$P
db$modules$annotations <- z
}

refs <- readRDS("results/inputs/reference_sets.rds")
db$gene_sets$sets <- c(refs$hallmark, refs$expanded, refs$cell)
