#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(optparse); library(Matrix)})

calculate_tmm <- function(counts, samples) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR is required for this optional sensitivity analysis; no installation is attempted automatically.")
  if (ncol(counts) != nrow(samples) || anyDuplicated(samples$sample_id) || anyNA(samples)) stop("Invalid pseudobulk sample metadata")
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
  if (any(!is.finite(values) | values < 0 | values != floor(values))) stop("TMM requires nonnegative integer raw counts")
  if (any(colSums(counts) != samples$source_library_size) || any(samples$source_library_size <= 0)) stop("Source library totals disagree with broad-gene counts")
  results <- list()
  for (cohort in unique(samples$cohort)) {
    idx <- which(samples$cohort == cohort)
    x <- as.matrix(counts[, idx, drop = FALSE])
    lib <- samples$source_library_size[idx]
    if (ncol(x) < 2) stop("At least two pseudobulks are required per TMM cohort")
    y <- edgeR::DGEList(counts = x, lib.size = lib)
    y <- edgeR::calcNormFactors(y, method = "TMM", logratioTrim = 0.3,
      sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10)
    factors <- y$samples$norm.factors
    if (any(!is.finite(factors) | factors <= 0) || abs(mean(log(factors))) > 1e-10) stop("Invalid TMM scaling factors")
    effective <- lib * factors
    # Check the denominator calculation against edgeR's public CPM implementation.
    stopifnot(isTRUE(all.equal(unname(edgeR::cpm(y, log = FALSE)),
      unname(sweep(x, 2, effective, "/") * 1e6), tolerance = 1e-12)))
    results[[cohort]] <- data.frame(sample_id = samples$sample_id[idx], cohort = cohort,
      source_library_size = lib, norm_factor = factors, effective_library_size = effective,
      genes_supplied = nrow(x), nonzero_genes_in_cohort = sum(rowSums(x) > 0))
  }
  do.call(rbind, results)
}

main <- function() {
  opts <- parse_args(OptionParser(option_list = list(
    make_option("--input-dir", type = "character", dest = "input_dir"),
    make_option("--outdir", type = "character")
  )))
  if (is.null(opts$input_dir) || is.null(opts$outdir)) stop("Supply --input-dir and a fresh --outdir")
  if (file.exists(opts$outdir) && (!dir.exists(opts$outdir) || length(list.files(opts$outdir, all.files = TRUE, no.. = TRUE)))) stop("Use an empty output directory")
  samples <- read.csv(file.path(opts$input_dir, "pseudobulk_samples.csv"), stringsAsFactors = FALSE)
  counts <- readMM(gzfile(file.path(opts$input_dir, "broad_counts.mtx.gz")))
  factors <- calculate_tmm(counts, samples)
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
  write.csv(factors, file.path(opts$outdir, "tmm_factors.csv"), row.names = FALSE)
  method <- list(package = "edgeR", version = as.character(packageVersion("edgeR")),
    method = "TMM", logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10,
    reference = "edgeR automatic reference per cohort", factor_centering = "geometric mean one per cohort",
    input = "Untransformed broad-source pseudobulk counts; no filterByExpr or HVG restriction",
    denominator = "source_library_size * norm_factor", cpm_check = "Matches edgeR::cpm(log=FALSE)",
    zero_handling = "edgeR removes all-zero genes and omits genes zero in either sample for pairwise TMM; standard TMM, not TMMwsp",
    limitation = "Descriptive sensitivity; no dispersion estimation or differential-expression fit; no absolute-transcript interpretation")
  jsonlite::write_json(method, file.path(opts$outdir, "tmm_method.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(opts$outdir, "sessionInfo.txt"))
  cat("TMM normalization completed for", length(unique(samples$cohort)), "cohorts\n")
}

if (sys.nframe() == 0L) main()
