#!/usr/bin/env Rscript

# Load function definitions only; do not execute the exporter or fit a model.
suppressPackageStartupMessages(library(CoGAPS))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "../.."))
definitions <- parse(file.path(root, "scripts/export_cogaps_diagnostics_r.R"))
for (expr in definitions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr)
}
md <- list(
  version = package_version("3.22.0"),
  params = CogapsParams(nPatterns = 6, nIterations = 2000, seed = 2, sparseOptimization = TRUE),
  chisq = as.numeric(40:1), atomsA = as.numeric(1:40), atomsP = as.numeric(41:80),
  totalUpdates = 123, totalRunningTime = 10, meanChiSq = 1
)
metrics <- list(
  status = "ok", language = "R", package = "CoGAPS", result_path = "fixture.rds",
  nPatterns = 6, nIterations = 2000, seed = 2, sparseOptimization = TRUE,
  outputFrequency = 100, totalUpdates = 123, totalRunningTime = 10, meanChiSq = 1,
  chisq_trace_length = 40, atomsA_trace_length = 40, atomsP_trace_length = 40
)
passed <- character()
run <- function(name, change = function(env) NULL, expected = "unresolved") {
  env <- new.env()
  env$md <- md
  env$metrics <- metrics
  change(env)
  warnings <- character()
  actual <- withCallingHandlers(build_trace_df(env$md, env$metrics, "fixture.rds"),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  stopifnot(identical(attr(actual, "trace_timing")$status, expected))
  if (expected == "unresolved") {
    stopifnot(length(warnings) == 1, all(is.na(actual$phase)), all(is.na(actual$iteration_in_phase)),
      all(is.na(actual$iteration_overall)), all(actual$phase_source == "unresolved"))
  } else {
    stopifnot(length(warnings) == 0, actual$phase[20] == "equilibration", actual$phase[21] == "sampling",
      actual$iteration_in_phase[1] == 100, actual$iteration_in_phase[20] == 2000,
      actual$iteration_in_phase[21] == 100, actual$iteration_overall[21] == env$metrics$nIterations + 100)
  }
  passed <<- c(passed, name)
}
run("verified complete run", expected = "reconstructed")
run("non-divisible iteration budget has no invented final status point", function(e) {
  e$md$params@nIterations <- 2050L; e$metrics$nIterations <- 2050
}, "reconstructed")
run("missing metrics", function(e) e$metrics <- NULL)
run("unknown saved version", function(e) e$md$version <- package_version("3.24.0"))
run("wrong result file", function(e) e$metrics$result_path <- "other.rds")
run("failed run", function(e) e$metrics$status <- "error")
run("wrong implementation", function(e) e$metrics$language <- "Python")
run("mismatched seed", function(e) e$metrics$seed <- 3)
run("mismatched pattern count", function(e) e$metrics$nPatterns <- 7)
run("mismatched iteration budget", function(e) e$metrics$nIterations <- 2100)
run("mismatched sparse setting", function(e) e$metrics$sparseOptimization <- FALSE)
run("mismatched total updates", function(e) e$metrics$totalUpdates <- 124)
run("mismatched runtime", function(e) e$metrics$totalRunningTime <- 11)
run("mismatched final statistic", function(e) e$metrics$meanChiSq <- 2)
run("missing output frequency", function(e) e$metrics$outputFrequency <- NULL)
run("disabled status output", function(e) e$metrics$outputFrequency <- 0)
run("fractional output frequency", function(e) e$metrics$outputFrequency <- 100.5)
run("non-numeric frequency", function(e) e$metrics$outputFrequency <- "100")
run("nonfinite frequency", function(e) e$metrics$outputFrequency <- Inf)
run("frequency beyond iteration budget", function(e) e$metrics$outputFrequency <- 3000)
run("wrong schedule despite forty rows", function(e) e$metrics$outputFrequency <- 200)
run("truncated trace", function(e) e$md$chisq <- e$md$chisq[-40])
run("unequal atom trace length", function(e) e$md$atomsP <- e$md$atomsP[-40])
run("mismatched recorded trace length", function(e) e$metrics$atomsA_trace_length <- 39)
run("nonfinite trace value", function(e) e$md$chisq[2] <- NA_real_)
run("fixed matrix mode", function(e) e$md$params@whichMatrixFixed <- "A")
stopifnot(nrow(build_trace_df(list())) == 0L)
passed <- c(passed, "empty trace is retained without invented rows")
cat(jsonlite::toJSON(list(tests_passed = length(passed), tests = passed), auto_unbox = TRUE, pretty = TRUE), "\n")
