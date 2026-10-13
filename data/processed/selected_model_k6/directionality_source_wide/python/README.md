# Pattern Gene Directionality Analysis

This folder contains a replicate-aware directionality check for the genes that define the chosen CoGAPS model patterns.

CoGAPS weights identify genes associated with each pattern, but they do not by themselves indicate whether those genes are up- or down-regulated. Direction was estimated separately from the underlying count data by comparing IFN-beta stimulated cells to control cells.

## Inputs

- processed AnnData: `preprocessed_cells_hvg3000.h5ad` (hash in `directionality_method.json`)
- chosen CoGAPS result: `cogaps_K6_seed2_iter2000.h5ad` (hash in `directionality_method.json`)
- top genes per pattern: `15`
- total unique genes tested: `60`

## Method

- extracted top CoGAPS genes from the chosen model A matrix
- aggregated raw counts from `.layers["counts"]` into pseudobulk profiles
- divided each gene count by the sum of verified source-wide totals for the same cells, then multiplied by one million; the denominator includes all 15,706 genes in the supplied source, not only the 3,000 modeled genes
- computed paired donor-level `log2(stim CPM + 0.5) - log2(ctrl CPM + 0.5)`
- reported eligible donors, positive/negative/zero counts, mean/median, agreement with the mean sign, and descriptive donor-support labels; no confidence category or significance test
- summarized direction globally across donors and within cell types where enough cells were available

## Main Result

- `Pattern5` is treated as the IFN-associated pattern for this run.
- `Pattern5`: 14 of 15 top genes have a positive mean stimulation contrast; consult the donor counts and sensitivity results before interpreting a direction.
- `Pattern3` has the next largest up-in-stimulation count among top pattern genes: 13 of 15.
- This supports using directionality as evidence for interpreting the IFN-associated pattern, while keeping non-IFN or secondary pattern labels more cautious.

## Caveats

- Candidate outputs only. These source-wide summaries do not by themselves assess normalization sensitivity; run the separate sensitivity/reporting step before lesson integration. No Wilcoxon/FDR or evidence-strength fields are carried forward.
- This is a targeted directionality check for pattern genes, not a genome-wide differential-expression analysis.
- The by-cell-type summaries are descriptive and depend on cell counts per donor/condition/cell-type stratum.
- The current dataset is an IFN-beta stimulation PBMC dataset; any dengue interpretation should be framed as biological relevance to IFN/inflammatory programs, not direct dengue infection evidence.

## Output Files

- `pattern_top_genes.csv`: top genes extracted from each CoGAPS pattern
- `pattern_gene_directionality_global.csv`: global paired pseudobulk direction by gene
- `pattern_gene_directionality_by_celltype.csv`: cell-type paired pseudobulk direction by gene
- `pattern_direction_summary.csv`: one-row summary per pattern
- `pattern_gene_directionality_pairs_global.csv`: donor-level global log2FC values
- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-level cell-type log2FC values
- `pattern_gene_directionality_heatmap.png`: heatmap of global top-gene directionality
