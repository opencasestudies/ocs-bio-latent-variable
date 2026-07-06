# Shared setup for the CoGAPS PBMC case study.
# Source this file before running later case-study chunks.

find_project_root_r <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(path, "index.qmd")) && dir.exists(file.path(path, "data"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the case-study project root.")
    }
    path <- parent
  }
}

resolve_first_existing_r <- function(...) {
  candidates <- c(...)
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  candidates[[1]]
}

ROOT_R <- find_project_root_r()
DATA_DIR_R <- file.path(ROOT_R, "data")
PROCESSED_DIR_R <- file.path(ROOT_R, "data", "processed")
INPUT_DIR_R <- file.path(PROCESSED_DIR_R, "input")
MODEL_SELECTION_DIR_R <- file.path(PROCESSED_DIR_R, "k_selection")
SELECTED_MODEL_DIR_R <- file.path(PROCESSED_DIR_R, "selected_model_k6")
FIGURES_DIR_R <- file.path(PROCESSED_DIR_R, "figures")
LEGACY_RESULTS_DIR_R <- file.path(ROOT_R, "data", "results")
SELECTED_RESULTS_DIR_R <- file.path(SELECTED_MODEL_DIR_R, "python")
PY_K6_RESULTS_DIR_R <- SELECTED_RESULTS_DIR_R
PY_RESULTS_DIR_R <- file.path(ROOT_R, "data", "results_python")
R_RESULTS_DIR_R <- file.path(ROOT_R, "data", "results_r")
R_K6_RESULTS_DIR_R <- file.path(SELECTED_MODEL_DIR_R, "r")
SELECTED_R_RESULTS_DIR_R <- R_K6_RESULTS_DIR_R
SELECTED_PYTHON_RESULTS_DIR_R <- SELECTED_RESULTS_DIR_R

CHOSEN_K_R <- 6L
CHOSEN_SEED_R <- 2L
CHOSEN_ITER_R <- 2000L
STIM_LABEL_R <- "stim"

ARTIFACT_MANIFEST_JSON_R <- file.path(DATA_DIR_R, "artifact_manifest.json")
K_SELECTION_DECISION_JSON_R <- file.path(MODEL_SELECTION_DIR_R, "k_selection_decision_manifest.json")
PREPROCESSED_H5AD_R <- file.path(INPUT_DIR_R, "preprocessed_cells_hvg3000.h5ad")
COGAPS_INPUT_H5AD_R <- file.path(INPUT_DIR_R, "cogaps_input_genesxcells_hvg3000_float64.h5ad")
PREPROCESS_CONFIG_JSON_R <- file.path(INPUT_DIR_R, "preprocess_config_hvg3000.json")

PHASE1_LOWERK_SUMMARY_CSV_R <- file.path(MODEL_SELECTION_DIR_R, "phase1_lowerk_summary.csv")
PHASE1B_LOWERK_LONG_SUMMARY_CSV_R <- file.path(MODEL_SELECTION_DIR_R, "phase1b_lowerk_long_summary.csv")
PHASE3_HIGHK_SUMMARY_CSV_R <- file.path(MODEL_SELECTION_DIR_R, "phase3_highk_2000_summary.csv")
CANDIDATE_PATTERN_SUMMARY_CSV_R <- file.path(MODEL_SELECTION_DIR_R, "candidate_pattern_summary.csv")
CANDIDATE_K7_ACTIVITY_MATCHES_CSV_R <- file.path(MODEL_SELECTION_DIR_R, "candidate_best_k7_activity_matches.csv")

R_CHOSEN_RESULT_RDS_R <- file.path(R_K6_RESULTS_DIR_R, "cogaps_K6_seed2_iter2000.rds")
R_CHOSEN_METRICS_JSON_R <- file.path(R_K6_RESULTS_DIR_R, "cogaps_K6_seed2_iter2000.metrics.json")
R_PATTERN_GENE_WEIGHTS_TSV_GZ_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_gene_weights.tsv.gz")
R_PATTERN_CELL_ACTIVITIES_TSV_GZ_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_cell_activities.tsv.gz")
R_PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_cell_activities_with_metadata.csv.gz")
R_PATTERN_TOP_GENES_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_top_genes.csv")
R_PATTERN_CORRELATIONS_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_correlations.csv")
R_PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_activity_by_celltype_condition.csv")
R_PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_activity_by_replicate_condition.csv")
R_PATTERN_SUMMARY_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_summary.csv")
R_DIRECTION_GLOBAL_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_gene_directionality_global.csv")
R_DIRECTION_CELLTYPE_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_gene_directionality_by_celltype.csv")
R_DIRECTION_SUMMARY_CSV_R <- file.path(R_K6_RESULTS_DIR_R, "pattern_direction_summary.csv")

PY_CHOSEN_RESULT_H5AD_R <- file.path(SELECTED_RESULTS_DIR_R, "cogaps_K6_seed2_iter2000.h5ad")
PY_CHOSEN_METRICS_JSON_R <- file.path(SELECTED_RESULTS_DIR_R, "cogaps_K6_seed2_iter2000.metrics.json")

HAS_R_CHOSEN_RESULT_R <- file.exists(R_CHOSEN_RESULT_RDS_R)
HAS_R_LIGHTWEIGHT_PATTERN_MATRICES_R <- file.exists(R_PATTERN_GENE_WEIGHTS_TSV_GZ_R) && file.exists(R_PATTERN_CELL_ACTIVITIES_TSV_GZ_R)
HAS_R_LIGHTWEIGHT_ACTIVITY_METADATA_R <- file.exists(R_PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ_R)
HAS_R_PATTERN_TOP_GENES_R <- file.exists(R_PATTERN_TOP_GENES_CSV_R)
HAS_R_PATTERN_CORRELATIONS_R <- file.exists(R_PATTERN_CORRELATIONS_CSV_R)
HAS_R_PATTERN_ACTIVITY_BY_CELLTYPE_R <- file.exists(R_PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV_R)
HAS_R_PATTERN_SUMMARY_R <- file.exists(R_PATTERN_SUMMARY_CSV_R)

R_ANALYSIS_MODE_R <- if (HAS_R_CHOSEN_RESULT_R) {
  "full_chosen_result"
} else if (HAS_R_LIGHTWEIGHT_PATTERN_MATRICES_R) {
  "lightweight_derived"
} else {
  "unavailable"
}

pattern_columns_r <- function(x) {
  pats <- grep("^Pattern_?[0-9]+$", x, value = TRUE)
  pats[order(as.integer(sub("^Pattern_?", "", pats)))]
}

canonicalize_pattern_names_r <- function(x) {
  sub("^Pattern_?([0-9]+)$", "Pattern\\1", x)
}

top_genes_from_matrix_r <- function(matrix, pattern, n = 15) {
  ordered <- sort(matrix[, pattern], decreasing = TRUE)
  ordered <- head(ordered, n)
  tibble::tibble(
    gene = names(ordered),
    weight = unname(ordered)
  )
}

load_lightweight_pattern_matrices_r <- function(weights_path, activities_path) {
  A <- readr::read_tsv(weights_path, show_col_types = FALSE) |>
    tibble::column_to_rownames("gene") |>
    as.data.frame(check.names = FALSE)
  P <- readr::read_tsv(activities_path, show_col_types = FALSE) |>
    tibble::column_to_rownames("cell_barcode") |>
    as.data.frame(check.names = FALSE)
  colnames(A) <- canonicalize_pattern_names_r(colnames(A))
  colnames(P) <- canonicalize_pattern_names_r(colnames(P))
  list(A = A, P = P)
}
