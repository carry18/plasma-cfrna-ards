options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})

out_root <- "results/figures"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

group_cols <- c(
  "Healthy controls" = "#8C9AA0",
  "Pneumonia without ARDS" = "#527D9A",
  "nonintubated ARDS" = "#93749E",
  "MV ARDS" = "#B38A56"
)
contrast_cols <- c(
  "ARDS overall vs pneumonia without ARDS" = "#2F5D7E",
  "nonintubated ARDS vs pneumonia without ARDS" = "#93749E",
  "MV ARDS vs pneumonia without ARDS" = "#B38A56",
  "MV vs nonintubated ARDS" = "#B36562"
)
signal_cols <- c(higher = "#B36562", lower = "#527D9A", neutral = "#A7AFB3")

theme_pub <- function(base_size = 7) theme_bw(base_family = "Arial", base_size = base_size) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#EDF0F2", linewidth = .2), panel.border = element_rect(colour = "#829099", linewidth = .35), legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_text(size = 6.5), plot.margin = margin(6, 6, 6, 6), strip.background = element_rect(fill = "#F4F6F7", colour = NA), strip.text = element_text(size = 7))
theme_set(theme_pub())

wrap_lab <- function(x, width = 30) vapply(x, function(z) paste(strwrap(z, width), collapse = "\n"), character(1))
clean_pathway <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("Il6", "IL6", x, fixed = TRUE)
  x <- gsub("Jak", "JAK", x, fixed = TRUE)
  x <- gsub("Stat3", "STAT3", x, fixed = TRUE)
  x <- gsub("Tnfa", "TNF-alpha", x, fixed = TRUE)
  x <- gsub("Nfkb", "NF-kB", x, fixed = TRUE)
  x <- gsub("E2f", "E2F", x, fixed = TRUE)
  x <- gsub("G2m", "G2M", x, fixed = TRUE)
  x <- gsub("Mtorc1", "mTORC1", x, fixed = TRUE)
  x <- gsub("Pi3k", "PI3K", x, fixed = TRUE)
  x
}

panelize <- function(p, tag) {
  p + labs(title = NULL, subtitle = NULL, caption = NULL, tag = tag) +
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 10),
      plot.tag.position = c(0, 1)
    )
}

save_figure <- function(p, number, height_mm = 180, width_mm = 183) {
  d <- file.path(out_root, paste0("Figure", number))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(d, paste0("Figure", number))
  ggsave(paste0(base, ".png"), p,
    width = width_mm, height = height_mm,
    units = "mm", dpi = 300, device = ragg::agg_png, bg = "white"
  )
  svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
  print(p)
  dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"),
    width = width_mm / 25.4,
    height = height_mm / 25.4, family = "Arial"
  )
  print(p)
  dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"),
    width = width_mm, height = height_mm,
    units = "mm", res = 600, compression = "lzw", background = "white"
  )
  print(p)
  dev.off()
  invisible(d)
}

write_source <- function(x, fig, name) {
  d <- file.path(out_root, paste0("Figure", fig), "source_data")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file.path(d, paste0(name, ".csv")), row.names = FALSE)
}

contrast_labels <- c(
  ARDS_equal_vs_pneumonia = "ARDS overall vs pneumonia without ARDS",
  nonMV_vs_pneumonia = "nonintubated ARDS vs pneumonia without ARDS",
  MV_vs_pneumonia = "MV ARDS vs pneumonia without ARDS",
  MV_vs_nonMV = "MV vs nonintubated ARDS"
)
