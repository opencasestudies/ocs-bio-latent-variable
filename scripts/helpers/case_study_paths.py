"""Shared setup for the CoGAPS PBMC case study.

Import this module before running later case-study chunks.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def find_project_root(start: Path) -> Path:
    for path in [start, *start.parents]:
        if (path / "index.qmd").exists() and (path / "data").exists():
            return path
    raise FileNotFoundError("Could not locate the case-study project root.")


def resolve_first_existing(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


ROOT = find_project_root(Path.cwd())
DATA_DIR = ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
INPUT_DIR = PROCESSED_DIR / "input"
MODEL_SELECTION_DIR = PROCESSED_DIR / "k_selection"
SELECTED_MODEL_DIR = PROCESSED_DIR / "selected_model_k6"
FIGURES_DIR = PROCESSED_DIR / "figures"
LEGACY_RESULTS_DIR = ROOT / "data" / "results"
SELECTED_RESULTS_DIR = SELECTED_MODEL_DIR / "python"
PY_K6_RESULTS_DIR = SELECTED_RESULTS_DIR
PY_RESULTS_DIR = ROOT / "data" / "results_python"
R_RESULTS_DIR = ROOT / "data" / "results_r"
R_K6_RESULTS_DIR = SELECTED_MODEL_DIR / "r"
SELECTED_R_RESULTS_DIR = R_K6_RESULTS_DIR
SELECTED_PYTHON_RESULTS_DIR = SELECTED_RESULTS_DIR

ARTIFACT_MANIFEST_JSON = DATA_DIR / "artifact_manifest.json"
K_SELECTION_DECISION_JSON = MODEL_SELECTION_DIR / "k_selection_decision_manifest.json"
PREPROCESSED_H5AD = INPUT_DIR / "preprocessed_cells_hvg3000.h5ad"
COGAPS_INPUT_H5AD = INPUT_DIR / "cogaps_input_genesxcells_hvg3000_float64.h5ad"
PREPROCESS_CONFIG_JSON = INPUT_DIR / "preprocess_config_hvg3000.json"

PY_CHOSEN_RESULT_H5AD = SELECTED_RESULTS_DIR / "cogaps_K6_seed2_iter2000.h5ad"
PY_CHOSEN_METRICS_JSON = SELECTED_RESULTS_DIR / "cogaps_K6_seed2_iter2000.metrics.json"
PY_SUMMARY_BY_K_CSV = resolve_first_existing(
    MODEL_SELECTION_DIR / "phase1b_lowerk_long_summary.csv",
    PY_RESULTS_DIR / "summary_by_K.csv",
    LEGACY_RESULTS_DIR / "summary_by_K.csv",
)
PHASE1_LOWERK_SUMMARY_CSV = MODEL_SELECTION_DIR / "phase1_lowerk_summary.csv"
PHASE1B_LOWERK_LONG_SUMMARY_CSV = MODEL_SELECTION_DIR / "phase1b_lowerk_long_summary.csv"
PHASE3_HIGHK_SUMMARY_CSV = MODEL_SELECTION_DIR / "phase3_highk_2000_summary.csv"
CANDIDATE_PATTERN_SUMMARY_CSV = MODEL_SELECTION_DIR / "candidate_pattern_summary.csv"
CANDIDATE_K7_ACTIVITY_MATCHES_CSV = MODEL_SELECTION_DIR / "candidate_best_k7_activity_matches.csv"
PY_PER_RUN_METRICS_CSV = resolve_first_existing(
    PY_RESULTS_DIR / "per_run_metrics.csv",
    LEGACY_RESULTS_DIR / "per_run_metrics.csv",
)
PY_DIRECTION_GLOBAL_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_gene_directionality_global.csv",
    PY_RESULTS_DIR / "pattern_gene_directionality_global.csv",
    LEGACY_RESULTS_DIR / "pattern_gene_directionality_global.csv",
)
PY_DIRECTION_CELLTYPE_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_gene_directionality_by_celltype.csv",
    PY_RESULTS_DIR / "pattern_gene_directionality_by_celltype.csv",
    LEGACY_RESULTS_DIR / "pattern_gene_directionality_by_celltype.csv",
)
PY_DIRECTION_SUMMARY_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_direction_summary.csv",
    PY_RESULTS_DIR / "pattern_direction_summary.csv",
    LEGACY_RESULTS_DIR / "pattern_direction_summary.csv",
)
PY_PATTERN_GENE_WEIGHTS_TSV_GZ = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_gene_weights.tsv.gz",
    PY_RESULTS_DIR / "pattern_gene_weights.tsv.gz",
    LEGACY_RESULTS_DIR / "pattern_gene_weights.tsv.gz",
)
PY_PATTERN_CELL_ACTIVITIES_TSV_GZ = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_cell_activities.tsv.gz",
    PY_RESULTS_DIR / "pattern_cell_activities.tsv.gz",
    LEGACY_RESULTS_DIR / "pattern_cell_activities.tsv.gz",
)
PY_PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_cell_activities_with_metadata.csv.gz",
    PY_RESULTS_DIR / "pattern_cell_activities_with_metadata.csv.gz",
    LEGACY_RESULTS_DIR / "pattern_cell_activities_with_metadata.csv.gz",
)
PY_PATTERN_TOP_GENES_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_top_genes.csv",
    PY_RESULTS_DIR / "pattern_top_genes.csv",
    LEGACY_RESULTS_DIR / "pattern_top_genes.csv",
)
PY_PATTERN_CORRELATIONS_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_correlations.csv",
    PY_RESULTS_DIR / "pattern_correlations.csv",
    LEGACY_RESULTS_DIR / "pattern_correlations.csv",
)
PY_PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_activity_by_celltype_condition.csv",
    PY_RESULTS_DIR / "pattern_activity_by_celltype_condition.csv",
    LEGACY_RESULTS_DIR / "pattern_activity_by_celltype_condition.csv",
)
PY_PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_activity_by_replicate_condition.csv",
    PY_RESULTS_DIR / "pattern_activity_by_replicate_condition.csv",
    LEGACY_RESULTS_DIR / "pattern_activity_by_replicate_condition.csv",
)
PY_PATTERN_SUMMARY_CSV = resolve_first_existing(
    SELECTED_RESULTS_DIR / "pattern_summary.csv",
    PY_RESULTS_DIR / "pattern_summary.csv",
    LEGACY_RESULTS_DIR / "pattern_summary.csv",
)

R_CHOSEN_RESULT_RDS = R_K6_RESULTS_DIR / "cogaps_K6_seed2_iter2000.rds"
R_CHOSEN_METRICS_JSON = R_K6_RESULTS_DIR / "cogaps_K6_seed2_iter2000.metrics.json"
R_PATTERN_GENE_WEIGHTS_TSV_GZ = R_K6_RESULTS_DIR / "pattern_gene_weights.tsv.gz"
R_PATTERN_CELL_ACTIVITIES_TSV_GZ = R_K6_RESULTS_DIR / "pattern_cell_activities.tsv.gz"
R_PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ = R_K6_RESULTS_DIR / "pattern_cell_activities_with_metadata.csv.gz"
R_PATTERN_TOP_GENES_CSV = R_K6_RESULTS_DIR / "pattern_top_genes.csv"
R_PATTERN_CORRELATIONS_CSV = R_K6_RESULTS_DIR / "pattern_correlations.csv"
R_PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV = R_K6_RESULTS_DIR / "pattern_activity_by_celltype_condition.csv"
R_PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV = R_K6_RESULTS_DIR / "pattern_activity_by_replicate_condition.csv"
R_PATTERN_SUMMARY_CSV = R_K6_RESULTS_DIR / "pattern_summary.csv"
R_DIRECTION_GLOBAL_CSV = R_K6_RESULTS_DIR / "pattern_gene_directionality_global.csv"
R_DIRECTION_CELLTYPE_CSV = R_K6_RESULTS_DIR / "pattern_gene_directionality_by_celltype.csv"
R_DIRECTION_SUMMARY_CSV = R_K6_RESULTS_DIR / "pattern_direction_summary.csv"

CHOSEN_K = 6
CHOSEN_SEED = 2
CHOSEN_ITER = 2000
STIM_LABEL = "stim"

HAS_CACHED_COGAPS_INPUT = COGAPS_INPUT_H5AD.exists()

HAS_PY_CHOSEN_RESULT = PY_CHOSEN_RESULT_H5AD.exists()
HAS_PY_LIGHTWEIGHT_PATTERN_MATRICES = (
    PY_PATTERN_GENE_WEIGHTS_TSV_GZ.exists() and PY_PATTERN_CELL_ACTIVITIES_TSV_GZ.exists()
)
HAS_PY_LIGHTWEIGHT_ACTIVITY_METADATA = PY_PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ.exists()
HAS_PY_PATTERN_TOP_GENES = PY_PATTERN_TOP_GENES_CSV.exists()
HAS_PY_PATTERN_CORRELATIONS = PY_PATTERN_CORRELATIONS_CSV.exists()
HAS_PY_PATTERN_ACTIVITY_BY_CELLTYPE = PY_PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV.exists()
HAS_PY_PATTERN_SUMMARY = PY_PATTERN_SUMMARY_CSV.exists()

HAS_R_CHOSEN_RESULT = R_CHOSEN_RESULT_RDS.exists()
HAS_R_LIGHTWEIGHT_PATTERN_MATRICES = (
    R_PATTERN_GENE_WEIGHTS_TSV_GZ.exists() and R_PATTERN_CELL_ACTIVITIES_TSV_GZ.exists()
)
HAS_R_LIGHTWEIGHT_ACTIVITY_METADATA = R_PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ.exists()
HAS_R_PATTERN_TOP_GENES = R_PATTERN_TOP_GENES_CSV.exists()
HAS_R_PATTERN_CORRELATIONS = R_PATTERN_CORRELATIONS_CSV.exists()
HAS_R_PATTERN_ACTIVITY_BY_CELLTYPE = R_PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV.exists()
HAS_R_PATTERN_SUMMARY = R_PATTERN_SUMMARY_CSV.exists()

PY_ANALYSIS_MODE = (
    "full_chosen_result"
    if HAS_PY_CHOSEN_RESULT
    else "lightweight_derived"
    if HAS_PY_LIGHTWEIGHT_PATTERN_MATRICES
    else "unavailable"
)
R_ANALYSIS_MODE = (
    "full_chosen_result"
    if HAS_R_CHOSEN_RESULT
    else "lightweight_derived"
    if HAS_R_LIGHTWEIGHT_PATTERN_MATRICES
    else "unavailable"
)

# Backward-compatible aliases used by the Python chunks below.
CHOSEN_RESULT_H5AD = PY_CHOSEN_RESULT_H5AD
CHOSEN_METRICS_JSON = PY_CHOSEN_METRICS_JSON
SUMMARY_BY_K_CSV = PY_SUMMARY_BY_K_CSV
PER_RUN_METRICS_CSV = PY_PER_RUN_METRICS_CSV
DIRECTION_GLOBAL_CSV = PY_DIRECTION_GLOBAL_CSV
DIRECTION_CELLTYPE_CSV = PY_DIRECTION_CELLTYPE_CSV
DIRECTION_SUMMARY_CSV = PY_DIRECTION_SUMMARY_CSV
PATTERN_GENE_WEIGHTS_TSV_GZ = PY_PATTERN_GENE_WEIGHTS_TSV_GZ
PATTERN_CELL_ACTIVITIES_TSV_GZ = PY_PATTERN_CELL_ACTIVITIES_TSV_GZ
PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ = PY_PATTERN_CELL_ACTIVITIES_WITH_METADATA_CSV_GZ
PATTERN_TOP_GENES_CSV = PY_PATTERN_TOP_GENES_CSV
PATTERN_CORRELATIONS_CSV = PY_PATTERN_CORRELATIONS_CSV
PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV = PY_PATTERN_ACTIVITY_BY_CELLTYPE_CONDITION_CSV
PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV = PY_PATTERN_ACTIVITY_BY_REPLICATE_CONDITION_CSV
PATTERN_SUMMARY_CSV = PY_PATTERN_SUMMARY_CSV
ANALYSIS_MODE = PY_ANALYSIS_MODE


def pattern_columns(df: pd.DataFrame) -> list[str]:
    return sorted(
        [c for c in df.columns if str(c).startswith("Pattern")],
        key=lambda value: int(str(value).replace("Pattern", "")),
    )


def top_genes_from_matrix(matrix: pd.DataFrame, pattern: str, n: int = 15) -> pd.DataFrame:
    top = matrix[pattern].sort_values(ascending=False).head(n)
    top_df = top.rename("weight").reset_index()
    top_df.columns = ["gene", "weight"]
    return top_df


def load_lightweight_pattern_matrices(weights_path: Path, activities_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    A = pd.read_csv(weights_path, sep="\t", compression="gzip").set_index("gene")
    P = pd.read_csv(activities_path, sep="\t", compression="gzip").set_index("cell_barcode")
    return A, P


def build_cogaps_input(preprocessed_cells: ad.AnnData) -> ad.AnnData:
    cogaps_input = preprocessed_cells.T.copy()
    X = cogaps_input.X.toarray() if sparse.issparse(cogaps_input.X) else cogaps_input.X
    cogaps_input.X = np.asarray(X, dtype=np.float64)
    cogaps_input.var = preprocessed_cells.obs.copy()
    cogaps_input.var_names = preprocessed_cells.obs_names.copy()
    cogaps_input.obs = preprocessed_cells.var.copy()
    cogaps_input.obs_names = preprocessed_cells.var_names.copy()
    return cogaps_input


__all__ = [
    "find_project_root",
    "resolve_first_existing",
    "pattern_columns",
    "top_genes_from_matrix",
    "load_lightweight_pattern_matrices",
    "build_cogaps_input",
    *[name for name in globals() if name.isupper()],
]
