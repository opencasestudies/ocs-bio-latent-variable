#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(CoGAPS)
  library(SummarizedExperiment)
  library(Matrix)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

lookup_count <- function(tbl, key) {
  if (key %in% names(tbl)) {
    as.integer(tbl[[key]])
  } else {
    0L
  }
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else "."
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
default_preprocessed <- file.path(ROOT, "data", "processed", "input", "preprocessed_cells_hvg3000.h5ad")
default_result_rds <- file.path(ROOT, "data", "processed", "selected_model_k6", "r", "cogaps_K6_seed2_iter2000.rds")
default_out_dir <- file.path(ROOT, "data", "processed", "selected_model_k6", "r")

option_list <- list(
  make_option("--preprocessed-h5ad", type = "character", dest = "preprocessed_h5ad", default = default_preprocessed, help = "Processed cells x genes AnnData (.h5ad) [default %default]"),
  make_option("--result-rds", type = "character", dest = "result_rds", default = default_result_rds, help = "Saved CoGAPS R result (.rds) [default %default]"),
  make_option("--normalization-manifest", type = "character", dest = "normalization_manifest", help = "Required manifest from scripts/directionality_normalization.py"),
  make_option("--outdir", type = "character", default = default_out_dir, help = "Output directory [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

PREPROCESSED_H5AD <- opt$preprocessed_h5ad
CHOSEN_RESULT_RDS <- opt$result_rds
OUT_DIR <- opt$outdir
NORMALIZATION_MANIFEST <- opt$normalization_manifest

TOP_N_GENES <- 15
MIN_CELLS_PER_PSEUDOBULK <- 20
CPM_PSEUDOCOUNT <- 0.5

pattern_columns <- function(x) {
  pats <- grep("^Pattern_?[0-9]+$", x, value = TRUE)
  pats[order(as.integer(sub("^Pattern_?", "", pats)))]
}

canonicalize_pattern_names <- function(x) {
  sub("^Pattern_?([0-9]+)$", "Pattern\\1", x)
}

ensure_pattern_matrix_orientation <- function(mat, expected_rows, expected_names) {
  if (nrow(mat) != expected_rows && ncol(mat) == expected_rows) {
    mat <- t(mat)
  }
  if (nrow(mat) != expected_rows) {
    stop(sprintf("Unexpected matrix shape. Expected %d rows, got %d.", expected_rows, nrow(mat)))
  }
  require_ids(rownames(mat), "model identifiers")
  if (!setequal(rownames(mat), expected_names)) {
    stop("Model identifiers differ from prepared identifiers; row order cannot substitute for names.")
  }
  mat[expected_names, , drop = FALSE]
}

summarize_donors <- function(values) {
  if (!is.numeric(values) || !is.null(dim(values)) || any(!is.finite(values))) stop("Donor contrasts must be a finite numeric vector")
  n <- length(values)
  positive <- sum(values > 0)
  negative <- sum(values < 0)
  zero <- sum(values == 0)
  avg <- if (n) mean(values) else NA_real_
  med <- if (n) median(values) else NA_real_
  support <- if (!n) "not estimable" else if (n == 1) "single donor pair; no across-donor replication" else if (positive > 0 && negative > 0) "mixed donor directions" else if (zero == n) "no observed contrast at this scale" else sprintf("%d of %d %s", max(positive, negative), n, if (positive) "positive" else "negative")
  agree <- if (!n || avg == 0) NA_real_ else if (avg > 0) positive else negative
  sign <- if (!n) "not estimable" else if (avg > 0) "positive" else if (avg < 0) "negative" else "zero"
  list(mean_log2fc_stim_vs_ctrl = avg, median_log2fc_stim_vs_ctrl = med,
    n_pairs = n, n_positive = positive, n_negative = negative, n_zero = zero,
    n_agree_mean = agree, agreement_with_mean = if (n) agree / n else NA_real_,
    mean_direction = sign, donor_support = support,
    direction = switch(sign, positive = "up_in_stim", negative = "down_in_stim", zero = "zero_mean", "not estimable" = "not_tested"),
    normalization_status = "not assessed; run normalization sensitivity")
}

extract_top_pattern_genes <- function(A, top_n = 15) {
  do.call(
    rbind,
    lapply(pattern_columns(colnames(A)), function(pattern) {
      ordered <- sort(A[, pattern], decreasing = TRUE)
      ordered <- head(ordered, top_n)
      data.frame(
        pattern = pattern,
        gene = names(ordered),
        rank = seq_along(ordered),
        cogaps_weight = unname(ordered),
        row.names = NULL
      )
    })
  )
}

get_counts_assay <- function(sce) {
  if (!"counts" %in% assayNames(sce)) {
    stop("Prepared counts assay is required; transformed X is not a raw-count fallback.")
  }
  counts <- assay(sce, "counts")
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
  if (any(!is.finite(values) | values < 0 | values != floor(values)) || sum(values) >= 2^53) {
    stop("Raw counts must be finite nonnegative integers within the exact-integer range.")
  }
  counts
}

require_ids <- function(ids, label) {
  if (is.null(ids) || !length(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop(paste(label, "must be unique and nonmissing"))
  }
}

file_record <- function(path) {
  list(filename = basename(path), size_bytes = unname(file.info(path)$size),
       sha256 = digest::digest(file = path, algo = "sha256"))
}

verify_file <- function(path, record) {
  actual <- file_record(path)
  if (!identical(actual$sha256, record$sha256) || actual$size_bytes != record$size_bytes) {
    stop(paste("Checksum/size mismatch:", path))
  }
}

load_normalization <- function(path, prepared_path, sce) {
  manifest <- jsonlite::read_json(path, simplifyVector = TRUE)
  if (!identical(manifest$schema_version, "source-wide-counts-v1") ||
      !isTRUE(manifest$retained_counts_match_source)) stop("Unsupported or unverified normalization manifest")
  verify_file(prepared_path, manifest$prepared)
  tables <- lapply(c("totals", "source_genes", "modeled_genes"), function(key) {
    table_path <- file.path(dirname(path), manifest[[key]]$filename)
    verify_file(table_path, manifest[[key]])
    read.csv(table_path, colClasses = "character", check.names = FALSE, na.strings = "")
  })
  names(tables) <- c("totals", "source_genes", "modeled_genes")
  require_ids(colnames(sce), "prepared cells")
  require_ids(rownames(sce), "prepared genes")
  require_ids(tables$totals$cell_id, "normalization cells")
  require_ids(tables$source_genes$gene, "normalization source genes")
  require_ids(tables$modeled_genes$gene, "normalization modeled genes")
  if (!setequal(tables$totals$cell_id, colnames(sce)) || ncol(sce) != manifest$n_cells ||
      !setequal(tables$modeled_genes$gene, rownames(sce)) ||
      !all(rownames(sce) %in% tables$source_genes$gene) ||
      nrow(sce) != manifest$modeled_gene_count || nrow(tables$source_genes) != manifest$source_gene_count) {
    stop("Normalization cell/gene universes disagree with prepared data")
  }
  totals <- tables$totals[match(colnames(sce), tables$totals$cell_id), ]
  hvg_totals <- as.numeric(colSums(get_counts_assay(sce)))
  if (!identical(as.numeric(totals$modeled_total_counts), hvg_totals)) {
    stop("Normalization modeled totals disagree with prepared counts")
  }
  for (column in c("replicate", "condition", "cell_type")) {
    values <- as.character(colData(sce)[[column]])
    if (length(values) != ncol(sce) || anyNA(values) || any(!nzchar(trimws(values))) ||
        !identical(totals[[column]], values)) stop(paste("Normalization metadata mismatch:", column))
  }
  if (!setequal(totals$condition, c("ctrl", "stim"))) stop("Expected ctrl and stim conditions")
  values <- as.numeric(totals$source_total_counts)
  if (any(!is.finite(values) | values <= 0 | values != floor(values) | values < hvg_totals) ||
      sum(values) >= 2^53) stop("Invalid source-wide library totals")
  list(totals = setNames(values, totals$cell_id), manifest = manifest)
}

group_eligibility <- function(obs, group_cols, min_cells) {
  groups <- aggregate(rep(1L, nrow(obs)), obs[group_cols], sum)
  names(groups)[ncol(groups)] <- "n_cells"
  groups$passes_min_cells <- groups$n_cells >= min_cells
  pair_cols <- setdiff(group_cols, "condition")
  keys <- do.call(paste, c(groups[pair_cols], sep = "\r"))
  groups$paired_eligible <- groups$passes_min_cells & ave(as.integer(groups$passes_min_cells), keys, FUN = sum) == 2
  groups$exclusion_reason <- ifelse(!groups$passes_min_cells, "below_min_cells",
    ifelse(groups$paired_eligible, "included", "counterpart_missing_or_below_min_cells"))
  groups
}

sum_gene_counts <- function(mat, gene_idx, cell_idx) {
  block <- mat[gene_idx, cell_idx, drop = FALSE]
  if (inherits(block, "sparseMatrix")) {
    return(as.numeric(Matrix::rowSums(block)))
  }
  rowSums(as.matrix(block))
}

group_library_size <- function(source_totals, cell_idx) {
  values <- source_totals[cell_idx]
  if (any(!is.finite(values) | values <= 0)) stop("Each source library total must be finite and positive")
  sum(values)
}

make_pseudobulk <- function(sce, genes, group_cols, min_cells, source_totals) {
  counts_mat <- get_counts_assay(sce)
  require_ids(colnames(sce), "prepared cells")
  require_ids(names(source_totals), "source-total cells")
  if (!setequal(names(source_totals), colnames(sce))) stop("Source totals must cover exactly the prepared cells")
  source_totals <- source_totals[colnames(sce)]
  gene_idx <- match(genes, rownames(sce))
  if (anyNA(gene_idx)) {
    stop("Some genes were not found in the processed SingleCellExperiment.")
  }

  obs <- as.data.frame(colData(sce))
  obs$.cell_index <- seq_len(ncol(sce))
  if (anyNA(obs[group_cols])) stop("Missing pseudobulk grouping metadata")

  keys <- do.call(paste, c(obs[group_cols], sep = "\r"))
  groups <- split(seq_len(nrow(obs)), keys)

  rows <- list()
  row_i <- 1L
  for (idx in groups) {
    group_df <- obs[idx, , drop = FALSE]
    if (nrow(group_df) < min_cells) {
      next
    }

    cell_idx <- group_df$.cell_index
    gene_counts <- sum_gene_counts(counts_mat, gene_idx, cell_idx)
    row <- as.list(group_df[1, group_cols, drop = FALSE])
    row$n_cells <- nrow(group_df)
    row$library_size_total_counts <- group_library_size(source_totals, cell_idx)
    row$library_size_pattern_genes <- sum(gene_counts)
    for (j in seq_along(genes)) {
      row[[genes[[j]]]] <- gene_counts[[j]]
    }
    rows[[row_i]] <- as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    row_i <- row_i + 1L
  }

  if (!length(rows)) {
    stop("No pseudobulk groups pass the minimum-cell rule")
  }
  do.call(rbind, rows)
}

add_cpm <- function(pseudobulk, genes) {
  out <- pseudobulk
  if (any(!is.finite(out$library_size_total_counts) | out$library_size_total_counts <= 0)) {
    stop("CPM requires finite positive source-wide library totals")
  }
  for (gene in genes) {
    count_col <- gene
    if (!count_col %in% names(out)) {
      stop(sprintf("Could not find pseudobulk count column for gene '%s'.", gene))
    }
    if (any(!is.finite(out[[count_col]]) | out[[count_col]] < 0 |
            out[[count_col]] > out$library_size_total_counts)) stop("Invalid pseudobulk counts")
    out[[paste0(gene, "__cpm")]] <- (out[[count_col]] / out$library_size_total_counts) * 1e6
  }
  out
}

paired_directionality <- function(pseudobulk_cpm, genes, pair_cols, context_cols, expected_strata = NULL) {
  if (anyNA(pseudobulk_cpm[c(pair_cols, "condition")])) stop("Missing donor or condition metadata")
  split_keys <- do.call(paste, c(pseudobulk_cpm[pair_cols], sep = "\r"))
  grouped <- split(seq_len(nrow(pseudobulk_cpm)), split_keys)

  pair_rows <- list()
  pair_i <- 1L

  for (idx in grouped) {
    group <- pseudobulk_cpm[idx, , drop = FALSE]
    if (anyDuplicated(group$condition) || !all(group$condition %in% c("ctrl", "stim"))) {
      stop("Expected at most one ctrl and one stim pseudobulk per donor stratum")
    }
    if (!all(c("ctrl", "stim") %in% group$condition)) {
      next
    }
    ctrl <- group[group$condition == "ctrl", , drop = FALSE][1, ]
    stim <- group[group$condition == "stim", , drop = FALSE][1, ]

    pair_info <- as.list(ctrl[1, pair_cols, drop = FALSE])
    for (gene in genes) {
      ctrl_cpm <- as.numeric(ctrl[[paste0(gene, "__cpm")]])
      stim_cpm <- as.numeric(stim[[paste0(gene, "__cpm")]])
      if (any(!is.finite(c(ctrl_cpm, stim_cpm))) || min(ctrl_cpm, stim_cpm) < 0) {
        stop("Paired CPM values must be finite and nonnegative")
      }
      log2fc <- log2(stim_cpm + CPM_PSEUDOCOUNT) - log2(ctrl_cpm + CPM_PSEUDOCOUNT)
      pair_rows[[pair_i]] <- c(
        context_cols,
        pair_info,
        list(
          gene = gene,
          ctrl_cpm = ctrl_cpm,
          stim_cpm = stim_cpm,
          log2fc_stim_vs_ctrl = log2fc,
          ctrl_n_cells = as.integer(ctrl$n_cells),
          stim_n_cells = as.integer(stim$n_cells)
        )
      )
      pair_i <- pair_i + 1L
    }
  }

  pair_df <- if (length(pair_rows)) do.call(rbind.data.frame, c(pair_rows, stringsAsFactors = FALSE)) else
    as.data.frame(setNames(replicate(length(c(names(context_cols), pair_cols, "gene", "ctrl_cpm", "stim_cpm", "log2fc_stim_vs_ctrl", "ctrl_n_cells", "stim_n_cells")), character(), simplify = FALSE),
      c(names(context_cols), pair_cols, "gene", "ctrl_cpm", "stim_cpm", "log2fc_stim_vs_ctrl", "ctrl_n_cells", "stim_n_cells")))
  numeric_cols <- c("ctrl_cpm", "stim_cpm", "log2fc_stim_vs_ctrl", "ctrl_n_cells", "stim_n_cells")
  for (col in numeric_cols) {
    pair_df[[col]] <- as.numeric(pair_df[[col]])
  }

  summary_group_cols <- unique(c("gene", setdiff(pair_cols, "replicate")))
  summary_rows <- list()
  summary_i <- 1L
  strata_cols <- setdiff(pair_cols, "replicate")
  strata <- if (is.null(expected_strata)) unique(pseudobulk_cpm[strata_cols]) else unique(expected_strata[strata_cols])
  contexts <- if (!length(strata_cols)) list(list()) else lapply(seq_len(nrow(strata)), function(i) as.list(strata[i, , drop = FALSE]))
  for (context in contexts) {
    for (gene in genes) {
      keep <- pair_df$gene == gene
      for (key in names(context)) keep <- keep & pair_df[[key]] == context[[key]]
      gene_df <- pair_df[keep, , drop = FALSE]
      if (anyDuplicated(gene_df$replicate)) stop("Duplicate donor contrasts")
      summary_rows[[summary_i]] <- c(context_cols, context, list(gene = gene), summarize_donors(as.numeric(gene_df$log2fc_stim_vs_ctrl)))
      summary_i <- summary_i + 1L
    }
  }

  summary_df <- do.call(rbind.data.frame, c(summary_rows, stringsAsFactors = FALSE))
  numeric_summary_cols <- c(
    "mean_log2fc_stim_vs_ctrl", "median_log2fc_stim_vs_ctrl", "n_pairs", "n_positive",
    "n_negative", "n_zero", "n_agree_mean", "agreement_with_mean"
  )
  for (col in numeric_summary_cols) {
    summary_df[[col]] <- as.numeric(summary_df[[col]])
  }

  list(pair_df = pair_df, summary_df = summary_df)
}

make_pattern_direction_summary <- function(global_summary) {
  do.call(
    rbind,
    lapply(split(global_summary, global_summary$pattern), function(df) {
      direction_counts <- table(df$direction)
      up_genes <- head(df$gene[df$direction == "up_in_stim"][order(df$mean_log2fc_stim_vs_ctrl[df$direction == "up_in_stim"], decreasing = TRUE)], 5)
      down_genes <- head(df$gene[df$direction == "down_in_stim"][order(df$mean_log2fc_stim_vs_ctrl[df$direction == "down_in_stim"], decreasing = FALSE)], 5)
      data.frame(
        pattern = df$pattern[[1]],
        n_top_genes = nrow(df),
        n_up_in_stim = lookup_count(direction_counts, "up_in_stim"),
        n_down_in_stim = lookup_count(direction_counts, "down_in_stim"),
        n_zero_mean = lookup_count(direction_counts, "zero_mean"),
        n_mixed_donor_directions = sum(df$donor_support == "mixed donor directions"),
        n_single_donor = sum(df$n_pairs == 1),
        n_not_estimable = sum(df$n_pairs == 0),
        median_mean_log2fc = median(df$mean_log2fc_stim_vs_ctrl),
        top_up_genes = paste(up_genes, collapse = ", "),
        top_down_genes = paste(down_genes, collapse = ", "),
        row.names = NULL
      )
    })
  )
}

make_pattern_heatmap <- function(global_summary, output_path) {
  ordered <- global_summary[order(global_summary$pattern, global_summary$rank), ]
  ordered$label <- paste0(ordered$pattern, ":", ordered$gene)
  p <- ggplot(ordered, aes(x = "stim vs ctrl\nmean paired log2FC", y = reorder(label, rev(seq_along(label))), fill = mean_log2fc_stim_vs_ctrl)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(x = NULL, y = NULL, title = "Directionality of top CoGAPS pattern genes", fill = "mean paired log2FC") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 6))
  ggsave(output_path, plot = p, width = 5, height = max(8, nrow(ordered) * 0.18), dpi = 200)
}

if (is.null(NORMALIZATION_MANIFEST)) {
  stop("Supply --normalization-manifest from scripts/directionality_normalization.py; no HVG-denominator fallback is permitted.")
}
if (file.exists(OUT_DIR) && (!dir.exists(OUT_DIR) || length(list.files(OUT_DIR, all.files = TRUE, no.. = TRUE)))) {
  stop("Use an empty output directory; existing files are protected.")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sce <- zellkonverter::readH5AD(PREPROCESSED_H5AD, reader = "R")
normalization <- load_normalization(NORMALIZATION_MANIFEST, PREPROCESSED_H5AD, sce)
result <- readRDS(CHOSEN_RESULT_RDS)
A <- as.matrix(getFeatureLoadings(result))
A <- ensure_pattern_matrix_orientation(A, nrow(sce), rownames(sce))
P <- ensure_pattern_matrix_orientation(as.matrix(getSampleFactors(result)), ncol(sce), colnames(sce))
if (any(!is.finite(A) | A < 0)) stop("Model gene weights must be finite and nonnegative")
if (is.null(colnames(A))) {
  colnames(A) <- sprintf("Pattern%d", seq_len(ncol(A)))
} else {
  colnames(A) <- canonicalize_pattern_names(colnames(A))
}

top_genes <- extract_top_pattern_genes(A, TOP_N_GENES)
genes <- unique(top_genes$gene)

global_pb <- add_cpm(
  make_pseudobulk(
    sce = sce,
    genes = genes,
    group_cols = c("replicate", "condition"),
    min_cells = MIN_CELLS_PER_PSEUDOBULK,
    source_totals = normalization$totals
  ),
  genes
)
global_direction <- paired_directionality(
  global_pb,
  genes = genes,
  pair_cols = c("replicate"),
  context_cols = list(contrast_scope = "global")
)
global_summary <- merge(top_genes, global_direction$summary_df, by = "gene", all.x = TRUE)
global_summary <- global_summary[order(global_summary$pattern, global_summary$rank), ]

by_cell_pb <- add_cpm(
  make_pseudobulk(
    sce = sce,
    genes = genes,
    group_cols = c("replicate", "condition", "cell_type"),
    min_cells = MIN_CELLS_PER_PSEUDOBULK,
    source_totals = normalization$totals
  ),
  genes
)
cell_direction <- paired_directionality(
  by_cell_pb,
  genes = genes,
  pair_cols = c("replicate", "cell_type"),
  context_cols = list(contrast_scope = "by_cell_type"),
  expected_strata = as.data.frame(colData(sce))["cell_type"]
)
cell_summary <- merge(top_genes, cell_direction$summary_df, by = "gene", all.x = TRUE)
cell_summary <- cell_summary[order(cell_summary$pattern, cell_summary$rank, cell_summary$cell_type), ]

pattern_summary <- make_pattern_direction_summary(global_summary)

write.csv(top_genes, file.path(OUT_DIR, "pattern_top_genes_directionality.csv"), row.names = FALSE)
write.csv(global_summary, file.path(OUT_DIR, "pattern_gene_directionality_global.csv"), row.names = FALSE)
write.csv(cell_summary, file.path(OUT_DIR, "pattern_gene_directionality_by_celltype.csv"), row.names = FALSE)
write.csv(pattern_summary, file.path(OUT_DIR, "pattern_direction_summary.csv"), row.names = FALSE)
write.csv(global_direction$pair_df, file.path(OUT_DIR, "pattern_gene_directionality_pairs_global.csv"), row.names = FALSE)
write.csv(cell_direction$pair_df, file.path(OUT_DIR, "pattern_gene_directionality_pairs_by_celltype.csv"), row.names = FALSE)

if (nrow(global_pb)) {
  write.csv(global_pb, file.path(OUT_DIR, "pseudobulk_counts_global.csv"), row.names = FALSE)
}
if (nrow(by_cell_pb)) {
  write.csv(by_cell_pb, file.path(OUT_DIR, "pseudobulk_counts_by_celltype.csv"), row.names = FALSE)
}

make_pattern_heatmap(global_summary, file.path(OUT_DIR, "pattern_gene_directionality_heatmap.png"))

for (scope in c("global", "by_celltype")) {
  cols <- c("replicate", "condition", if (scope == "by_celltype") "cell_type")
  write.csv(group_eligibility(as.data.frame(colData(sce)), cols, MIN_CELLS_PER_PSEUDOBULK),
    file.path(OUT_DIR, paste0("pseudobulk_eligibility_", scope, ".csv")), row.names = FALSE)
}
method <- list(
  normalization = normalization$manifest,
  model = file_record(CHOSEN_RESULT_RDS), script = file_record(script_path),
  numerator = "Sum of prepared counts for one gene in the exact grouped cells",
  denominator = "Sum of source-wide per-cell raw totals for those same cells",
  grouping_global = c("replicate", "condition"),
  grouping_cell_type = c("replicate", "condition", "cell_type"),
  min_cells_per_pseudobulk = MIN_CELLS_PER_PSEUDOBULK,
  pseudocount = CPM_PSEUDOCOUNT, pseudocount_unit = "CPM",
  contrast = "log2(stim_CPM + 0.5) - log2(ctrl_CPM + 0.5)",
  donor_weighting = "Equal weight among eligible donor-paired contrasts",
  status = "Candidate outputs; not approved for replacing lesson results",
  reporting_version = "donor-support-v1",
  caution = "Descriptive donor counts and mean-sign agreement only; no evidence-strength category or Wilcoxon/FDR. Source-only exports mark sensitivity not assessed until the sensitivity step.",
  outputs = lapply(list.files(OUT_DIR, full.names = TRUE), file_record)
)
jsonlite::write_json(method, file.path(OUT_DIR, "directionality_method.json"), auto_unbox = TRUE, pretty = TRUE)

cat("Wrote R directionality outputs to:", OUT_DIR, "\n")
