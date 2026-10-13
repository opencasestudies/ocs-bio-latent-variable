#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "../.."))
source(file.path(root, "scripts/calculate_tmm_factors.R"))
counts <- cbind(c(4, 6, 20, 30), c(8, 12, 40, 60), c(16, 24, 80, 120))
samples <- data.frame(sample_id = c("a", "b", "c"), cohort = "fixture", source_library_size = colSums(counts))
result <- calculate_tmm(counts, samples)
stopifnot(max(abs(result$norm_factor - 1)) < 1e-12)
stopifnot(max(abs(result$effective_library_size - colSums(counts))) < 1e-10)
permutation <- c(3, 1, 2)
other <- calculate_tmm(counts[, permutation], samples[permutation, ])
stopifnot(max(abs(other$norm_factor[match(samples$sample_id, other$sample_id)] - result$norm_factor)) < 1e-12)
reject <- function(expr) stopifnot(tryCatch({force(expr); FALSE}, error = function(e) TRUE))
for (value in c(-1, 0.5, NA_real_, Inf)) {
  bad <- counts
  bad[1, 1] <- value
  reject(calculate_tmm(bad, samples))
}
bad <- samples
bad$source_library_size[1] <- 1
reject(calculate_tmm(counts, bad))
bad <- samples
bad$sample_id[1] <- "b"
reject(calculate_tmm(counts, bad))
reject(calculate_tmm(counts[, 1, drop = FALSE], samples[1, ]))
cat("10 TMM checks passed: proportional libraries, effective sizes, permutation, four invalid counts, wrong totals, duplicate samples, one sample.\n")
