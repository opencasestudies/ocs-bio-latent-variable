# Case Study Data

This folder contains the GitHub-safe artifacts used by the Bio-OCS CoGAPS case study, plus documentation for optional full-local reproduction files.

The main learner path does not rerun the full model sweep. Instead, it loads a processed teaching input, precomputed model-selection summaries, and lightweight exports from the selected `K = 6`, `seed = 2`, `n_iter = 2000`, sparse-on CoGAPS model.

## Manifest

- `artifact_manifest.json`: machine-readable inventory of the files used by the case study, including file sizes, checksums, and optional external/local artifacts.

## `processed/input/`

- `preprocessed_cells_hvg3000.h5ad`: GitHub-safe processed cells x genes AnnData object used for the main learner workflow.
- `preprocess_config_hvg3000.json`: preprocessing configuration used to create the cached object.
- `cogaps_input_genesxcells_hvg3000_float64.h5ad`: optional local-only dense genes x cells matrix for full Python/PyCoGAPS reruns. This file is ignored by Git and should be supplied through a local download or external artifact host if needed.

## `processed/k_selection/`

These files summarize the precomputed sweep evidence used for rank selection. Learners inspect these summaries rather than rerunning the HPC sweep.

- `k_selection_decision_manifest.json`: compact statement of the current selected model and evidence files.
- `phase1_lowerk_summary.csv`: lower-rank sweep summary.
- `phase1b_lowerk_long_summary.csv`: longer-iteration lower-rank summary.
- `phase3_highk_2000_summary.csv`: high-rank sweep summary.
- `candidate_pattern_summary.csv`: pattern-level summaries for candidate models.
- `candidate_best_k7_activity_matches.csv`: compact K6-to-K7 reference matching summary.

## `processed/selected_model_k6/r/`

Selected R/Bioconductor CoGAPS `K = 6` sparse-on artifacts used by the primary case-study path.

GitHub-safe files include:

- `cogaps_K6_seed2_iter2000.metrics.json`: selected R model metadata.
- `cogaps_K6_seed2_iter2000.diagnostics.json`: selected R diagnostic summary.
- `cogaps_K6_seed2_iter2000.trace.csv`: chi-square and atom-count trace.
- `cogaps_K6_seed2_iter2000.snapshot_summary.csv`: snapshot stability summary.
- `cogaps_K6_seed2_iter2000.pattern_uncertainty_summary.csv`: pattern uncertainty summary.
- `cogaps_K6_seed2_iter2000.pump_stat.tsv.gz`: pump-sampling output.
- `cogaps_K6_seed2_iter2000.mean_pattern_assignment.tsv.gz`: mean pattern-assignment output.
- `pattern_gene_weights.tsv.gz`: lightweight exported gene-weight (`A`) matrix.
- `pattern_gene_weight_sd.tsv.gz`: gene-weight standard deviations.
- `pattern_cell_activities.tsv.gz`: lightweight exported cell-activity (`P`) matrix.
- `pattern_cell_activity_sd.tsv.gz`: cell-activity standard deviations.
- `pattern_cell_activities_with_metadata.csv.gz`: cell-level pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern.
- `pattern_correlations.csv`: correlations between pattern activity and stimulation.
- `pattern_activity_by_celltype_condition.csv`: compact activity summaries by cell type and condition.
- `pattern_activity_by_replicate_condition.csv`: donor-aware activity summaries.
- `pattern_summary.csv`: compact per-pattern interpretation table.
- `pattern_gene_directionality_global.csv`: targeted global directionality for top pattern genes.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for top pattern genes.
- `pattern_direction_summary.csv`: one-row directionality summary per pattern.
- `pattern_top_genes_directionality.csv`: narrower top-gene subset used in directionality summaries.
- `pattern_gene_directionality_pairs_global.csv`: donor-paired global log2FC table.
- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-paired cell-type log2FC table.
- `pattern_gene_directionality_heatmap.png`: heatmap summarizing top pattern-gene directionality.
- `pseudobulk_counts_global.csv`: global pseudobulk counts used for directionality summaries.
- `pseudobulk_counts_by_celltype.csv`: cell-type pseudobulk counts used for directionality summaries.

Optional full-local files, ignored by Git:

- `cogaps_K6_seed2_iter2000.rds`: full selected R CoGAPS object.
- `cogaps_K6_seed2_iter2000.checkpoint.out`: checkpoint from the selected heavy diagnostic run.

## `processed/selected_model_k6/python/`

Selected Python/PyCoGAPS-compatible `K = 6` sparse-on artifacts used in the parallel Python tabs.

The Python files mirror the R artifact families where possible: metrics, pattern matrices, activity summaries, top genes, and directionality summaries. The optional full `.h5ad` model object is ignored by Git and is only needed for full local inspection or regeneration of Python-derived exports.

## `processed/figures/`

Learner-facing figure assets copied from the selected K6 analysis. These support Quarto rendering and keep figures near the data products that generated them.

## `external/`

This directory documents optional large files that may be supplied through Zenodo, another external artifact host, or local private storage. The ordinary GitHub-rendered case study should not require these files.

## Legacy Reference Outputs

Older K7, distributed, runtime-comparison, and developer-dashboard outputs may still exist elsewhere under `data/` for provenance. They are not part of the selected learner-path artifact contract unless a section explicitly labels them as reference or sensitivity evidence.
