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

## Trace Timing Provenance

The trace has 40 saved status points, not 40 individual sampler iterations. The saved R object does not attach phase or iteration labels to those points. Its parameters record 2,000 iterations per phase; the matching run metrics record output every 100 iterations. Inspection of CoGAPS 3.22.0 source at commit `4118fd66c028954ddddce9455c0aaa25e2e58968` verifies that status values are appended after each output-frequency multiple, first during equilibration and then during sampling, with the iteration counter reset between phases.

The four added columns explicitly identify the reconstructed timing:

| Column | Meaning |
| --- | --- |
| `phase` | `equilibration` or `sampling`, after the run-record checks pass. |
| `iteration_in_phase` | Completed iteration within that phase: 100, 200, ..., 2,000. |
| `iteration_overall` | Completed iteration across both phases: 100, 200, ..., 4,000. |
| `phase_source` | `reconstructed_from_run_metrics`, not directly recorded phase labels. |

Points 1-20 belong to equilibration; points 21-40 belong to sampling. Point 21 is sampling iteration 100, or iteration 2,100 overall. `relative_progress` is preserved for compatibility but is only the saved point index divided by the point count; it is not an independently recorded iteration.

The diagnostics JSON records the reconstruction method, source revision, consistency checks, and MD5 file identities for the saved result and run metrics. The original chi-square/atom-count values, uncertainty summaries, snapshots, and model are unchanged. The final 20 chi-square values range from 222,414,208 to 222,438,496; this does not establish convergence.

`scripts/export_cogaps_diagnostics_r.R` automatically looks for the `.metrics.json` beside the chosen `.rds`; `--run-metrics` can specify another location. Reconstruction is limited to the verified CoGAPS 3.22.0 convention with matching run identifiers/statistics and a complete output schedule. If the record is missing, inconsistent, or unsupported, the raw trace is still exported, but phase/iteration fields remain missing, `phase_source` is `unresolved`, and a warning explains why. The exporter never labels an arbitrary trace by simply dividing its rows in half.

The selected run's console log was not located. The reconstruction is grounded in its saved object, matching metrics, recorded phase-specific snapshots, and the matching source implementation, not a claimed log-based recovery. The learner plots now consume these columns with checks for resolved and consistent timing. The full-run view uses overall iterations; the sampling-only view uses within-phase iterations and an explicitly expanded vertical scale. Both retain the original chi-square values.
