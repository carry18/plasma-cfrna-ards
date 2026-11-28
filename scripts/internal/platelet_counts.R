# Recorded platelet-count comparisons; one value per participant.
summarise_platelet_counts <- function(clinical) {
  stopifnot(all(c("Sample", "group", "platelet") %in% names(clinical)),
            !anyDuplicated(clinical$Sample), is.numeric(clinical$platelet))
  groups <- lapply(c("A", "B", "C"), function(group) {
    clinical$platelet[clinical$group == group & is.finite(clinical$platelet)]
  })
  names(groups) <- c("A", "B", "C")
  groups$AC <- c(groups$A, groups$C)
  comparisons <- list(A_B = c("A", "B"), A_C = c("A", "C"),
                      B_C = c("B", "C"), B_AC = c("B", "AC"))
  result <- do.call(rbind, lapply(names(comparisons), function(comparison) {
    pair <- comparisons[[comparison]]
    left <- groups[[pair[1]]]
    right <- groups[[pair[2]]]
    stopifnot(length(left) > 0L, length(right) > 0L)
    ql <- quantile(left, c(.25, .5, .75), names = FALSE, type = 7)
    qr <- quantile(right, c(.25, .5, .75), names = FALSE, type = 7)
    data.frame(
      comparison = comparison, n_left = length(left), n_right = length(right),
      median_left_10e9_L = ql[2], Q1_left = ql[1], Q3_left = ql[3],
      median_right_10e9_L = qr[2], Q1_right = qr[1], Q3_right = qr[3],
      nominal_P = wilcox.test(left, right, alternative = "two.sided",
                             exact = FALSE, correct = TRUE)$p.value
    )
  }))
  result$Holm_P_four_comparisons <- p.adjust(result$nominal_P, method = "holm")
  rownames(result) <- NULL
  result
}
