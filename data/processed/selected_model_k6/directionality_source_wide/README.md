# Source-Wide Directionality Results

This bundle contains corrected expression contrasts for the selected K6 models,
explicit donor-support summaries, and three-method normalization comparisons.
The core lesson reads these corrected reports, and the Visualization source
tables and gene-contrast figures use the same estimates and donor reporting.
The Full Reproduction Guide provides the matching saved-figure, source-count,
TMM-comparison and optional full-model export commands using runtime 0.3.2.
Superseded exports have been removed from the model-summary folders. Keep this
bundle's folder structure: its table schemas and reporting rules differ from
the retired exports.

## Files

| Folder | Contents |
| --- | --- |
| `normalization/` | Verified source-wide cell totals, source and modeled gene identifiers, and input fingerprints. |
| `r/` and `python/` | Source-wide stimulation exports, paired donor values, eligibility records, method metadata, and their saved heatmaps. These source-only summaries mark normalization sensitivity as unassessed. |
| `sensitivity/` | Three-method donor contrasts, summaries, comparison flags, eligible-donor coverage, and language-specific reports. |

For learner tables, use `sensitivity/r_gene_report_global.csv`,
`r_gene_report_by_celltype.csv`, and `r_gene_report_identity.csv`, or the files
with the `python_` prefix. These reports join source-wide estimates to the
normalization comparison. `eligible_donors_by_scope.csv` gives the donor coverage;
`normalization_comparison.csv` retains all three method means.

Source-wide CPM uses the total raw counts across the supplied object's 15,706
genes, not only the 3,000 modeled genes. The log2 contrasts use a 0.5-CPM
pseudocount and equal weight among eligible paired donors. Identity contrasts
pool all non-focal cells within donor and condition, then average the two
condition contrasts within each donor before summarizing donors. The model
inputs and fitted gene weights/cell activities have not changed.

The primary reports include means, medians, eligible pairs, positive/negative/zero
counts, agreement with the mean sign, donor descriptions, and a separate
normalization flag. Same mean sign across methods is not a significance test,
agreement among donors, or evidence of absolute transcript changes. A single
donor pair does not provide across-donor replication. TMM is a sensitivity
comparison, not a preferred answer selected after inspecting the results.

## Integrity And Import Check

From the repository root, this read-only check uses the Python standard library:

```bash
python3 scripts/reproduction/check_directionality_bundle.py
```

It checks checksums, byte sizes, CSV schemas, reporting fields, component
consistency, and donor/normalization labels. It does not require source H5ADs,
full models, edgeR, or a model fit. Add
`--prepared data/processed/input/preprocessed_cells_hvg3000.h5ad` to verify the
included prepared file's identity as well.

The bundle manifest records all calculation payloads and their schemas. The
original generation manifests and outputs are byte-preserved, including their
candidate-status notes; those notes describe their generation stage rather than
publication approval. The top-level `integration_state` records the local
lesson and reproduction integration; it does not certify a public release.
The repository-wide inventory also includes this bundle and identifies optional
archive downloads. Check its listed files without downloading optional models:

```bash
python3 scripts/reproduction/check_artifact_manifest.py
```
