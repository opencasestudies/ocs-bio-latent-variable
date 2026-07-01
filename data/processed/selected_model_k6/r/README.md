# Selected R K6 CoGAPS Artifacts

This folder contains the primary R/Bioconductor CoGAPS artifacts for the selected `K = 6`, `seed = 2`, `n_iter = 2000`, sparse-on model.

The ordinary case-study path uses the lightweight files in this directory. The optional full `.rds` model object and checkpoint are ignored by Git and are only needed for full local inspection, derivative regeneration, or methods debugging.

Pattern labels are run-specific. In this selected local R run, the IFN-associated pattern is `Pattern4`; prose should use "the IFN-associated pattern" unless it is discussing this exact saved run.

## Lightweight Files

- `cogaps_K6_seed2_iter2000.metrics.json`: selected R model metadata.
- `cogaps_K6_seed2_iter2000.diagnostics.json`: compact diagnostic summary.
- `cogaps_K6_seed2_iter2000.trace.csv`: chi-square and atom-count trace.
- `cogaps_K6_seed2_iter2000.snapshot_summary.csv`: snapshot stability summary.
- `cogaps_K6_seed2_iter2000.pattern_uncertainty_summary.csv`: pattern uncertainty summary.
- `pattern_gene_weights.tsv.gz`: exported gene-weight (`A`) matrix.
- `pattern_cell_activities.tsv.gz`: exported cell-activity (`P`) matrix.
- `pattern_cell_activities_with_metadata.csv.gz`: cell-level pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern.
- `pattern_correlations.csv`: pattern correlations with stimulation.
- `pattern_activity_by_celltype_condition.csv`: activity summaries by cell type and condition.
- `pattern_activity_by_replicate_condition.csv`: donor-aware activity summaries.
- `pattern_summary.csv`: compact per-pattern interpretation table.
- `pattern_gene_directionality_global.csv`: targeted global directionality for top pattern genes.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for top pattern genes.
- `pattern_direction_summary.csv`: one-row directionality summary per pattern.

## Optional Full-Local Files

- `cogaps_K6_seed2_iter2000.rds`: full selected R CoGAPS object.
- `cogaps_K6_seed2_iter2000.checkpoint.out`: checkpoint from the selected heavy diagnostic run.

These files are not required for the main learner path.
