# 输入：03网络及对应输入。
# 输出：350轮重采样检查点。
# 在项目根目录运行；前置步骤：01、02、03。
suppressPackageStartupMessages(library(WGCNA))
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
out <- "results/analysis/stability"
dir.create(file.path(out, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
new <- readRDS("results/analysis/wgcna/analysis_checkpoint_before_external.rds")
a <- list(md = new$metadata, bc = new$discovery, e = new$expression, dat = t(new$expression[new$network_genes, new$discovery, drop = FALSE]), proj = new$scores, module_sets = new$module_sets, power = new$power)
md <- a$md[a$bc, ]
stopifnot(identical(md$Sample, rownames(a$dat)))
strata <- split(seq_len(nrow(md)), interaction(md$group, md$Batch, drop = TRUE))
all_dat <- t(a$e[, a$bc])
ref <- a$proj[md$Sample, , drop = FALSE]
groups <- as.character(md$group)
delta <- function(s, idx) mean(s[idx][groups[idx] == "B"]) - mean(s[idx][groups[idx] == "C"])
score_pc <- function(genes, idx) {
  x <- all_dat[, genes, drop = FALSE]
  mu <- colMeans(x[idx, , drop = FALSE])
  ss <- apply(x[idx, , drop = FALSE], 2, sd)
  if (any(!is.finite(ss) | ss <= 0)) stop("Zero/non-finite variance in module")
  z <- sweep(sweep(x, 2, mu, "-"), 2, ss, "/")
  pc <- prcomp(z[idx, , drop = FALSE], center = FALSE, scale. = FALSE, rank. = 1)
  w <- pc$rotation[, 1]
  if (stats::cor(pc$x[, 1], rowMeans(z[idx, , drop = FALSE])) < 0) w <- -w
  s <- as.numeric(z %*% w)
  s / sd(s[idx])
}
run_one <- function(mode, fraction, iteration) {
  # Distinct seeds across random-draw scenarios; leave-one-out indices are deterministic.
  seed <- 20261017L + iteration + ifelse(fraction == .9, 1000L, 0L) + ifelse(mode == "reselected", 2000L, 0L)
  set.seed(seed)
  idx <- sort(unlist(lapply(strata, function(v) sample(v, max(1, round(length(v) * fraction)), replace = FALSE))))
  if (mode == "loo") idx <- setdiff(seq_len(nrow(md)), iteration)
  omitted <- setdiff(seq_len(nrow(md)), idx)
  genes <- colnames(a$dat)
  if (mode == "reselected") {
    mg <- apply(all_dat[idx, , drop = FALSE], 2, mad)
    genes <- names(sort(mg, decreasing = TRUE))[seq_len(min(5000, sum(mg > 0)))]
  }
  warnings <- character()
  start <- proc.time()[3]
  ans <- tryCatch(withCallingHandlers(
    {
      net <- blockwiseModules(all_dat[idx, genes, drop = FALSE],
        power = a$power,
        networkType = "signed", TOMType = "signed", minModuleSize = 40, mergeCutHeight = .25,
        numericLabels = TRUE, pamRespectsDendro = FALSE, maxBlockSize = 6000,
        randomSeed = seed, verbose = 0, nThreads = 1
      )
      labs <- setNames(net$colors, genes)
      mods <- sort(unique(labs[labs > 0]))
      matches <- list()
      overlaps <- list()
      memberships <- list()
      for (nm in names(a$module_sets)) {
        original <- a$module_sets[[nm]]
        ov <- do.call(rbind, lapply(mods, function(k) {
          gg <- names(labs)[labs == k]
          hit <- length(intersect(original, gg))
          data.frame(original = nm, label = k, size = length(gg), intersection = hit, jaccard = hit / length(union(original, gg)), recall = hit / length(original), precision = hit / length(gg))
        }))
        if (is.null(ov)) stop("No non-grey modules recovered")
        best <- ov[order(-ov$jaccard, ov$label), ][1, ]
        gg <- names(labs)[labs == best$label]
        sc <- score_pc(gg, idx)
        best$score_cor_retained <- stats::cor(sc[idx], ref[idx, nm])
        best$score_cor_omitted <- if (length(omitted) > 1) stats::cor(sc[omitted], ref[omitted, nm]) else NA_real_
        best$effect_retained <- delta(sc, idx)
        best$effect_omitted <- delta(sc, omitted)
        best$effect_reference <- delta(ref[, nm], seq_len(nrow(md)))
        best$direction_retained <- sign(best$effect_retained) == sign(best$effect_reference)
        best$gene_universe_recall <- mean(original %in% genes)
        matches[[nm]] <- best
        overlaps[[nm]] <- ov
        memberships[[nm]] <- gg
      }
      list(status = "success", matches = do.call(rbind, matches), overlaps = do.call(rbind, overlaps), memberships = memberships, labels = labs)
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ), error = function(e) list(status = "failure", error = conditionMessage(e)))
  c(list(
    mode = mode, fraction = fraction, iteration = iteration, seed = seed,
    input_fingerprint = input_fingerprint, retained = md$Sample[idx], omitted = md$Sample[omitted], warnings = warnings, seconds = unname(proc.time()[3] - start)
  ), ans)
}
input_fingerprint <- digest::digest(list(metadata = md, expression = all_dat, reference_scores = ref, network_genes = colnames(a$dat), module_sets = a$module_sets, power = a$power, algorithm = body(run_one)), algo = "sha256")
limit <- 100L
scenarios <- data.frame(mode = c("fixed", "fixed", "reselected", "loo"), fraction = c(.8, .9, .8, 1))
for (j in seq_len(nrow(scenarios))) {
  for (i in seq_len(limit)) {
    mode <- scenarios$mode[j]
    frac <- scenarios$fraction[j]
    if (mode == "loo" && i > nrow(md)) next
    f <- file.path(out, "checkpoints", sprintf("%s_%02d_%03d.rds", mode, round(frac * 100), i))
    if (file.exists(f)) {
      saved <- readRDS(f)
      if (identical(saved$status, "success") && identical(saved$input_fingerprint, input_fingerprint)) next
      stop("Stale or failed checkpoint; use a clean output directory: ", f)
    }
    # One atomic directory lock per draw prevents duplicate work across runners.
    lock <- paste0(f, ".lock")
    if (!dir.create(lock, showWarnings = FALSE)) next
    if (file.exists(f)) {
      unlink(lock, recursive = TRUE)
      next
    }
    z <- run_one(mode, frac, i)
    saveRDS(z, f)
    unlink(lock, recursive = TRUE)
    cat(sprintf("%s %.0f%% %03d %s %.1fs\n", mode, frac * 100, i, z$status, z$seconds))
    flush.console()
    gc(FALSE)
  }
}
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

files <- list.files(file.path(out, "checkpoints"), pattern = "[.]rds$", full.names = TRUE)
stopifnot(length(files) == 350L, all(vapply(files, function(f) identical(readRDS(f)$status, "success"), logical(1))))
