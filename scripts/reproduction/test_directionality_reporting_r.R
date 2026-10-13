#!/usr/bin/env Rscript

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", arg[[1]])), "../.."))
for (expr in parse(file.path(root, "scripts/pattern_directionality_analysis_r.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr)
}
CPM_PSEUDOCOUNT <- 0.5
fixtures <- list(none = numeric(), one = 8, one_zero = 0, all_zero = c(0, 0),
  opposite_positive = c(10, -1, -1), opposite_negative = c(-10, 1, 1),
  positive_zero = c(2, 2, 0), negative_zero = c(-2, 0), tied_zero_mean = c(-1, -1, 1, 1),
  tiny_positive = c(1e-12, 1e-12))
expected <- c("not estimable", "single donor pair; no across-donor replication", "single donor pair; no across-donor replication",
  "no observed contrast at this scale", "mixed donor directions", "mixed donor directions", "2 of 3 positive", "1 of 2 negative",
  "mixed donor directions", "2 of 2 positive")
results <- lapply(fixtures, summarize_donors)
stopifnot(identical(unname(vapply(results, `[[`, character(1), "donor_support")), expected))
stopifnot(is.na(results$none$agreement_with_mean), is.na(results$one_zero$agreement_with_mean),
  is.na(results$all_zero$n_agree_mean), is.na(results$tied_zero_mean$agreement_with_mean),
  results$opposite_positive$agreement_with_mean == 1/3, results$opposite_negative$agreement_with_mean == 1/3,
  results$positive_zero$agreement_with_mean == 2/3)
for (bad in c(NA_real_, Inf, -Inf)) {
  stopifnot(inherits(try(summarize_donors(c(1, bad)), silent = TRUE), "try-error"))
}
pb <- data.frame(replicate = "d1", condition = "ctrl", cell_type = "A", n_cells = 20, g__cpm = 1)
result <- paired_directionality(pb, "g", c("replicate", "cell_type"), list())
stopifnot(nrow(result$pair_df) == 0, result$summary_df$donor_support == "not estimable")
pb <- rbind(pb, transform(pb, condition = "stim"))
result <- paired_directionality(pb, "g", c("replicate", "cell_type"), list(), data.frame(cell_type = c("A", "absent")))
stopifnot(result$summary_df$donor_support[result$summary_df$cell_type == "absent"] == "not estimable")
stopifnot(inherits(try(paired_directionality(rbind(pb, pb[1, ]), "g", c("replicate", "cell_type"), list()), silent = TRUE), "try-error"))
cat(jsonlite::toJSON(list(fixture_cases = length(fixtures), additional_checks = 12, results = results), pretty = TRUE, auto_unbox = TRUE), "\n")
