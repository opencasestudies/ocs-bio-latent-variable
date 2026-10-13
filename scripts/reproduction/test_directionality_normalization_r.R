#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "../.."))
# Load actual production functions without running the exporter.
for (expr in parse(file.path(root, "scripts/pattern_directionality_analysis_r.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr)
}
CPM_PSEUDOCOUNT <- 0.5
STRONG_ABS_LOG2FC <- 0.5
MODERATE_ABS_LOG2FC <- 0.25
HIGH_CONSISTENCY <- 0.75
passed <- character()
reject <- function(label, expr) {
  result <- tryCatch({force(expr); FALSE}, error = function(e) TRUE)
  stopifnot(result)
  passed <<- c(passed, label)
}
counts <- Matrix(matrix(c(1, 2, 2, 5, 2, 1, 4, 3), nrow = 2), sparse = TRUE)
dimnames(counts) <- list(c("g-1", "g2"), letters[1:4])
sce <- SummarizedExperiment(assays = list(counts = counts, X = counts / 2),
  colData = DataFrame(replicate = rep("d1", 4), condition = c("ctrl", "ctrl", "stim", "stim"), cell_type = rep("B", 4)))
totals <- setNames(c(6, 14, 13, 27), letters[1:4])
pb <- add_cpm(make_pseudobulk(sce, c("g-1", "g2"), c("replicate", "condition"), 2, rev(totals)), c("g-1", "g2"))
stopifnot(identical(pb$library_size_total_counts, c(20, 40)), identical(pb[["g-1"]], c(3, 6)),
  identical(pb[["g-1__cpm"]], c(150000, 150000)), "g-1" %in% names(pb))
pairs <- paired_directionality(pb, "g-1", "replicate", list())
stopifnot(pairs$pair_df$log2fc_stim_vs_ctrl == 0)
passed <- c(passed, "hand-calculated source CPM, hyphenated gene names, and reordered barcodes")
missing <- sce
assay(missing, "counts") <- NULL
reject("no transformed-X fallback", get_counts_assay(missing))
for (value in c(-1, 0.5, NA_real_, Inf, 2^53)) {
  bad <- sce
  assay(bad, "counts")[1, 1] <- value
  reject(paste("invalid count", value), get_counts_assay(bad))
}
for (ids in list(c("a", "a"), c("a", ""), c("a", NA_character_), NULL)) {
  reject("invalid identifiers", require_ids(ids, "fixture"))
}
for (value in c(0, -1, NA_real_, Inf)) {
  bad <- pb
  bad$library_size_total_counts[1] <- value
  reject(paste("invalid denominator", value), add_cpm(bad, "g-1"))
}
reject("mismatched source-total cells", make_pseudobulk(sce, "g-1", c("replicate", "condition"), 2, totals[-1]))
reject("no eligible groups", make_pseudobulk(sce, "g-1", c("replicate", "condition"), 3, totals))
reject("duplicate donor condition", paired_directionality(rbind(pb, pb[1, ]), "g-1", "replicate", list()))
bad <- pb
bad[["g-1__cpm"]][1] <- NA_real_
reject("nonfinite paired CPM", paired_directionality(bad, "g-1", "replicate", list()))
bad <- sce
colData(bad)$condition[1] <- NA_character_
reject("missing grouping metadata", make_pseudobulk(bad, "g-1", c("replicate", "condition"), 2, totals))
mat <- matrix(1:4, 2, dimnames = list(c("b", "a"), c("Pattern1", "Pattern2")))
stopifnot(identical(rownames(ensure_pattern_matrix_orientation(mat, 2, c("a", "b"))), c("a", "b")))
rownames(mat) <- NULL
reject("unnamed model rows are not guessed", ensure_pattern_matrix_orientation(mat, 2, c("a", "b")))
obs <- as.data.frame(colData(sce))
obs$cell_type[4] <- "T"
coverage <- group_eligibility(obs, c("replicate", "condition", "cell_type"), 2)
stopifnot(!any(coverage$paired_eligible), sum(coverage$exclusion_reason == "below_min_cells") == 2,
  sum(coverage$exclusion_reason == "counterpart_missing_or_below_min_cells") == 1)
passed <- c(passed, "explicit pair eligibility and exclusion reasons")
cat(jsonlite::toJSON(list(tests_passed = length(passed), tests = passed), pretty = TRUE, auto_unbox = TRUE), "\n")
