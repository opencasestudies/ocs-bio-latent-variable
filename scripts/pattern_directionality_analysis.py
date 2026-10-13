#!/usr/bin/env python3
"""Assess stim-vs-ctrl directionality for Python/PyCoGAPS pattern genes.

This script answers a reviewer-facing question that CoGAPS alone does not:
the CoGAPS A matrix identifies genes that define each pattern, but expression
direction needs to be estimated from the underlying count data.

The selected Python `.h5ad` is an optional full-model input. The Full Reproduction
Guide sets COGAPS_OUT_DIR to a fresh `learner_outputs/recalculated_expression/python/`
folder for this script's tables and heatmap; existing outputs are protected.
These exports use source-wide CPM but do not yet assess normalization sensitivity.
The main learner path instead uses included, sensitivity-qualified reports from
`data/processed/selected_model_k6/directionality_source_wide/sensitivity/` and does
not require the full model object.

The analysis uses replicate-aware pseudobulk summaries:

- extract top genes from the chosen CoGAPS result
- aggregate raw counts by replicate + condition, globally and within cell type
- library-size normalize to CPM with verified source-wide totals
- compute paired log2(stim CPM + pseudocount) - log2(ctrl CPM + pseudocount)
- summarize direction and consistency across donors
"""

from __future__ import annotations

from pathlib import Path
import os

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from directionality_normalization import (
    fresh_output, group_eligibility, load_normalization, prepared_counts,
    require_ids, write_method_manifest,
)
from directionality_reporting import stimulation_summary


ROOT = Path(__file__).resolve().parents[1]
PREPROCESSED_H5AD = Path(
    os.environ.get(
        "COGAPS_PREPROCESSED_H5AD",
        ROOT / "data/processed/input/preprocessed_cells_hvg3000.h5ad",
    )
)
CHOSEN_RESULT_H5AD = Path(
    os.environ.get(
        "COGAPS_RESULT_H5AD",
        ROOT / "data/processed/selected_model_k6/python/cogaps_K6_seed2_iter2000.h5ad",
    )
)
OUT_DIR = Path(os.environ.get("COGAPS_OUT_DIR", ROOT / "data/processed/selected_model_k6/python"))
IFN_PATTERN = os.environ.get("COGAPS_IFN_PATTERN", "").strip()
NORMALIZATION_MANIFEST = os.environ.get("COGAPS_NORMALIZATION_MANIFEST", "").strip()

TOP_N_GENES = 15
MIN_CELLS_PER_PSEUDOBULK = 20
CPM_PSEUDOCOUNT = 0.5


def pattern_columns(df: pd.DataFrame) -> list[str]:
    return sorted(
        [c for c in df.columns if str(c).startswith("Pattern")],
        key=lambda s: int(str(s).replace("Pattern", "")),
    )


def get_counts_matrix(adata: ad.AnnData):
    return prepared_counts(adata)


def subset_counts(X, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    block = X[rows][:, cols]
    if sparse.issparse(block):
        return np.asarray(block.sum(axis=0)).ravel()
    return np.asarray(block).sum(axis=0)


def group_library_size(source_totals, rows: np.ndarray) -> float:
    values = np.asarray(source_totals, dtype=np.float64)[rows]
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("Each source library total must be finite and positive")
    return float(values.sum())


def extract_top_pattern_genes(result: ad.AnnData, top_n: int) -> pd.DataFrame:
    rows = []
    patterns = pattern_columns(result.obs)
    for pattern in patterns:
        weights = result.obs[pattern].sort_values(ascending=False).head(top_n)
        for rank, (gene, weight) in enumerate(weights.items(), start=1):
            rows.append(
                {
                    "pattern": pattern,
                    "gene": str(gene),
                    "rank": rank,
                    "cogaps_weight": float(weight),
                }
            )
    return pd.DataFrame(rows)


def make_pseudobulk(
    adata: ad.AnnData,
    genes: list[str],
    group_cols: list[str],
    min_cells: int,
    source_totals: pd.Series,
) -> pd.DataFrame:
    missing = sorted(set(genes) - set(adata.var_names.astype(str)))
    if missing:
        raise ValueError(f"Pattern genes missing from processed AnnData: {missing}")

    X = get_counts_matrix(adata)
    cells = require_ids(adata.obs_names, "prepared cells")
    ids = require_ids(source_totals.index, "source-total cells")
    if set(ids) != set(cells):
        raise ValueError("Source totals must cover exactly the prepared cells")
    source_totals = source_totals.loc[cells]
    gene_indices = np.asarray([adata.var_names.get_loc(gene) for gene in genes])
    obs = adata.obs.copy()
    obs["_row_index"] = np.arange(adata.n_obs)
    if obs[group_cols].isna().any().any():
        raise ValueError("Missing pseudobulk grouping metadata")

    rows = []
    for group_values, group_df in obs.groupby(group_cols, observed=True):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        n_cells = len(group_df)
        if n_cells < min_cells:
            continue

        cell_indices = group_df["_row_index"].to_numpy()
        gene_counts = subset_counts(X, cell_indices, gene_indices)
        library_size = group_library_size(source_totals, cell_indices)

        row = dict(zip(group_cols, group_values))
        row["n_cells"] = int(n_cells)
        row["library_size_total_counts"] = library_size
        row["library_size_pattern_genes"] = float(gene_counts.sum())
        for gene, count in zip(genes, gene_counts):
            row[gene] = float(count)
        rows.append(row)

    if not rows:
        raise ValueError("No pseudobulk groups pass the minimum-cell rule")
    return pd.DataFrame(rows)


def add_cpm(pseudobulk: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    out = pseudobulk.copy()
    count_matrix = out[genes].to_numpy(dtype=float)
    library_sizes = out["library_size_total_counts"].to_numpy(dtype=float)
    if not np.isfinite(library_sizes).all() or (library_sizes <= 0).any():
        raise ValueError("CPM requires finite positive source-wide library totals")
    if (not np.isfinite(count_matrix).all() or (count_matrix < 0).any()
            or (count_matrix > library_sizes[:, None]).any()):
        raise ValueError("Invalid pseudobulk counts")
    cpm = count_matrix / library_sizes[:, None] * 1_000_000
    for i, gene in enumerate(genes):
        out[f"{gene}__cpm"] = cpm[:, i]
    return out


def paired_directionality(
    pseudobulk_cpm: pd.DataFrame,
    genes: list[str],
    pair_cols: list[str],
    context_cols: dict[str, str],
    expected_strata: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-pair log2FC rows and one summary row per gene."""
    pair_rows = []
    summary_rows = []

    group_cols = [*pair_cols]
    if pseudobulk_cpm[group_cols + ["condition"]].isna().any().any():
        raise ValueError("Missing donor or condition metadata")
    grouped = pseudobulk_cpm.groupby(group_cols, observed=True, dropna=False)
    for pair_values, group in grouped:
        if not isinstance(pair_values, tuple):
            pair_values = (pair_values,)
        by_condition = group.set_index("condition", drop=False)
        if by_condition.index.has_duplicates or not set(by_condition.index).issubset({"ctrl", "stim"}):
            raise ValueError("Expected at most one ctrl and one stim pseudobulk per donor stratum")
        if "ctrl" not in by_condition.index or "stim" not in by_condition.index:
            continue
        pair_info = dict(zip(pair_cols, pair_values))
        ctrl = by_condition.loc["ctrl"]
        stim = by_condition.loc["stim"]
        for gene in genes:
            ctrl_cpm = float(ctrl[f"{gene}__cpm"])
            stim_cpm = float(stim[f"{gene}__cpm"])
            if not np.isfinite([ctrl_cpm, stim_cpm]).all() or min(ctrl_cpm, stim_cpm) < 0:
                raise ValueError("Paired CPM values must be finite and nonnegative")
            log2fc = np.log2(stim_cpm + CPM_PSEUDOCOUNT) - np.log2(ctrl_cpm + CPM_PSEUDOCOUNT)
            pair_rows.append(
                {
                    **context_cols,
                    **pair_info,
                    "gene": gene,
                    "ctrl_cpm": ctrl_cpm,
                    "stim_cpm": stim_cpm,
                    "log2fc_stim_vs_ctrl": float(log2fc),
                    "ctrl_n_cells": int(ctrl["n_cells"]),
                    "stim_n_cells": int(stim["n_cells"]),
                }
            )

    pair_df = pd.DataFrame(pair_rows, columns=[*context_cols, *pair_cols, "gene", "ctrl_cpm", "stim_cpm",
                                              "log2fc_stim_vs_ctrl", "ctrl_n_cells", "stim_n_cells"])
    summary_group_cols = ["gene", *[c for c in pair_cols if c != "replicate"]]
    strata_cols = summary_group_cols[1:]
    strata = (expected_strata if expected_strata is not None else pseudobulk_cpm[strata_cols]).drop_duplicates()
    contexts = strata.to_dict("records") if strata_cols else [{}]
    for context in contexts:
        for gene in genes:
            mask = pair_df.gene.eq(gene)
            for key, value in context.items():
                mask &= pair_df[key].eq(value)
            group = pair_df.loc[mask]
            if group.replicate.duplicated().any():
                raise ValueError("Duplicate donor contrasts")
            summary_rows.append({**context_cols, **context, "gene": gene,
                                 **stimulation_summary(group.log2fc_stim_vs_ctrl.to_numpy(float))})
    return pair_df, pd.DataFrame(summary_rows)


def merge_pattern_annotations(summary: pd.DataFrame, top_genes: pd.DataFrame) -> pd.DataFrame:
    merged = top_genes.merge(summary, on="gene", how="left")
    sort_cols = [c for c in ["pattern", "rank", "cell_type"] if c in merged.columns]
    return merged.sort_values(sort_cols).reset_index(drop=True)


def make_pattern_heatmap(global_summary: pd.DataFrame, output_path: Path) -> None:
    ordered = global_summary.sort_values(["pattern", "rank"])
    labels = [f"{row.pattern}:{row.gene}" for row in ordered.itertuples()]
    values = ordered["mean_log2fc_stim_vs_ctrl"].to_numpy(dtype=float)[:, None]

    height = max(8, len(labels) * 0.18)
    fig, ax = plt.subplots(figsize=(5, height))
    vmax = max(1.0, float(np.nanmax(np.abs(values))))
    im = ax.imshow(values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks([0])
    ax.set_xticklabels(["stim vs ctrl\nmean paired log2FC"])
    ax.set_title("Directionality of top CoGAPS pattern genes")
    fig.colorbar(im, ax=ax, fraction=0.08, pad=0.04, label="mean paired log2FC")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def make_pattern_direction_summary(global_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pattern, df in global_summary.groupby("pattern", observed=True):
        direction_counts = df["direction"].value_counts().to_dict()
        up_genes = (
            df.loc[df["direction"] == "up_in_stim"]
            .sort_values("mean_log2fc_stim_vs_ctrl", ascending=False)
            .head(5)["gene"]
        )
        down_genes = (
            df.loc[df["direction"] == "down_in_stim"]
            .sort_values("mean_log2fc_stim_vs_ctrl", ascending=True)
            .head(5)["gene"]
        )
        rows.append(
            {
                "pattern": pattern,
                "n_top_genes": int(len(df)),
                "n_up_in_stim": int(direction_counts.get("up_in_stim", 0)),
                "n_down_in_stim": int(direction_counts.get("down_in_stim", 0)),
                "n_zero_mean": int(direction_counts.get("zero_mean", 0)),
                "n_mixed_donor_directions": int(df.donor_support.eq("mixed donor directions").sum()),
                "n_single_donor": int(df.n_pairs.eq(1).sum()),
                "n_not_estimable": int(df.n_pairs.eq(0).sum()),
                "median_mean_log2fc": float(df["mean_log2fc_stim_vs_ctrl"].median()),
                "top_up_genes": ", ".join(up_genes),
                "top_down_genes": ", ".join(down_genes),
            }
        )
    return pd.DataFrame(rows).sort_values("pattern")


def write_readme(
    output_path: Path,
    n_genes: int,
    global_summary: pd.DataFrame,
    pattern_summary: pd.DataFrame,
) -> None:
    if IFN_PATTERN and IFN_PATTERN in set(pattern_summary["pattern"]):
        ifn_pattern = IFN_PATTERN
    else:
        ifn_pattern = (
            pattern_summary
            .sort_values(["n_up_in_stim", "median_mean_log2fc"], ascending=False)
            .iloc[0]["pattern"]
        )
    ifn_row = pattern_summary.loc[pattern_summary["pattern"] == ifn_pattern].iloc[0]
    secondary = (
        pattern_summary
        .loc[pattern_summary["pattern"] != ifn_pattern]
        .sort_values(["n_up_in_stim", "median_mean_log2fc"], ascending=False)
    )
    secondary_row = secondary.iloc[0] if not secondary.empty else None

    lines = [
        "# Pattern Gene Directionality Analysis",
        "",
        "This folder contains a replicate-aware directionality check for the genes that define the chosen CoGAPS model patterns.",
        "",
        "CoGAPS weights identify genes associated with each pattern, but they do not by themselves indicate whether those genes are up- or down-regulated. Direction was estimated separately from the underlying count data by comparing IFN-beta stimulated cells to control cells.",
        "",
        "## Inputs",
        "",
        f"- processed AnnData: `{PREPROCESSED_H5AD.name}` (hash in `directionality_method.json`)",
        f"- chosen CoGAPS result: `{CHOSEN_RESULT_H5AD.name}` (hash in `directionality_method.json`)",
        f"- top genes per pattern: `{TOP_N_GENES}`",
        f"- total unique genes tested: `{n_genes}`",
        "",
        "## Method",
        "",
        "- extracted top CoGAPS genes from the chosen model A matrix",
        "- aggregated raw counts from `.layers[\"counts\"]` into pseudobulk profiles",
        "- divided each gene count by the sum of verified source-wide totals for the same cells, then multiplied by one million; the denominator includes all 15,706 genes in the supplied source, not only the 3,000 modeled genes",
        "- computed paired donor-level `log2(stim CPM + 0.5) - log2(ctrl CPM + 0.5)`",
        "- reported eligible donors, positive/negative/zero counts, mean/median, agreement with the mean sign, and descriptive donor-support labels; no confidence category or significance test",
        "- summarized direction globally across donors and within cell types where enough cells were available",
        "",
        "## Main Result",
        "",
        f"- `{ifn_pattern}` is treated as the IFN-associated pattern for this run.",
        f"- `{ifn_pattern}`: {int(ifn_row.n_up_in_stim)} of {int(ifn_row.n_top_genes)} top genes have a positive mean stimulation contrast; consult the donor counts and sensitivity results before interpreting a direction.",
        (
            f"- `{secondary_row.pattern}` has the next largest up-in-stimulation count among top pattern genes: "
            f"{int(secondary_row.n_up_in_stim)} of {int(secondary_row.n_top_genes)}."
            if secondary_row is not None else
            "- No secondary pattern summary is available."
        ),
        "- This supports using directionality as evidence for interpreting the IFN-associated pattern, while keeping non-IFN or secondary pattern labels more cautious.",
        "",
        "## Caveats",
        "",
        "- Candidate outputs only. These source-wide summaries do not by themselves assess normalization sensitivity; run the separate sensitivity/reporting step before lesson integration. No Wilcoxon/FDR or evidence-strength fields are carried forward.",
        "- This is a targeted directionality check for pattern genes, not a genome-wide differential-expression analysis.",
        "- The by-cell-type summaries are descriptive and depend on cell counts per donor/condition/cell-type stratum.",
        "- The current dataset is an IFN-beta stimulation PBMC dataset; any dengue interpretation should be framed as biological relevance to IFN/inflammatory programs, not direct dengue infection evidence.",
        "",
        "## Output Files",
        "",
        "- `pattern_top_genes.csv`: top genes extracted from each CoGAPS pattern",
        "- `pattern_gene_directionality_global.csv`: global paired pseudobulk direction by gene",
        "- `pattern_gene_directionality_by_celltype.csv`: cell-type paired pseudobulk direction by gene",
        "- `pattern_direction_summary.csv`: one-row summary per pattern",
        "- `pattern_gene_directionality_pairs_global.csv`: donor-level global log2FC values",
        "- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-level cell-type log2FC values",
        "- `pattern_gene_directionality_heatmap.png`: heatmap of global top-gene directionality",
        "",
    ]
    output_path.write_text("\n".join(lines))


def main() -> None:
    if not NORMALIZATION_MANIFEST:
        raise ValueError("Set COGAPS_NORMALIZATION_MANIFEST to the manifest produced by scripts/directionality_normalization.py; no HVG-denominator fallback is permitted")
    fresh_output(OUT_DIR)

    adata = ad.read_h5ad(PREPROCESSED_H5AD)
    source_totals, normalization = load_normalization(NORMALIZATION_MANIFEST, PREPROCESSED_H5AD, adata)
    result = ad.read_h5ad(CHOSEN_RESULT_H5AD)
    if (set(require_ids(result.obs_names, "model genes")) != set(adata.var_names)
            or set(require_ids(result.var_names, "model cells")) != set(adata.obs_names)):
        raise ValueError("Model identifiers differ from the prepared gene/cell sets")
    weights = result.obs[pattern_columns(result.obs)].to_numpy(dtype=float)
    if weights.shape[1] == 0 or not np.isfinite(weights).all() or (weights < 0).any():
        raise ValueError("Model gene weights must be finite and nonnegative")

    top_genes = extract_top_pattern_genes(result, TOP_N_GENES)
    genes = sorted(top_genes["gene"].unique())

    global_pb = add_cpm(
        make_pseudobulk(
            adata=adata,
            genes=genes,
            group_cols=["replicate", "condition"],
            min_cells=MIN_CELLS_PER_PSEUDOBULK,
            source_totals=source_totals,
        ),
        genes,
    )
    global_pairs, global_gene_summary = paired_directionality(
        global_pb,
        genes,
        pair_cols=["replicate"],
        context_cols={"contrast_scope": "global"},
    )
    global_summary = merge_pattern_annotations(global_gene_summary, top_genes)

    by_cell_pb = add_cpm(
        make_pseudobulk(
            adata=adata,
            genes=genes,
            group_cols=["replicate", "condition", "cell_type"],
            min_cells=MIN_CELLS_PER_PSEUDOBULK,
            source_totals=source_totals,
        ),
        genes,
    )
    cell_pairs, cell_gene_summary = paired_directionality(
        by_cell_pb,
        genes,
        pair_cols=["replicate", "cell_type"],
        context_cols={"contrast_scope": "by_cell_type"},
        expected_strata=adata.obs[["cell_type"]],
    )
    cell_summary = (
        top_genes.merge(cell_gene_summary, on="gene", how="left")
        .sort_values(["pattern", "rank", "cell_type"])
        .reset_index(drop=True)
    )

    pattern_summary = make_pattern_direction_summary(global_summary)

    top_genes.to_csv(OUT_DIR / "pattern_top_genes.csv", index=False)
    global_summary.to_csv(OUT_DIR / "pattern_gene_directionality_global.csv", index=False)
    cell_summary.to_csv(OUT_DIR / "pattern_gene_directionality_by_celltype.csv", index=False)
    pattern_summary.to_csv(OUT_DIR / "pattern_direction_summary.csv", index=False)
    global_pairs.to_csv(OUT_DIR / "pattern_gene_directionality_pairs_global.csv", index=False)
    cell_pairs.to_csv(OUT_DIR / "pattern_gene_directionality_pairs_by_celltype.csv", index=False)
    global_pb.to_csv(OUT_DIR / "pseudobulk_counts_global.csv", index=False)
    by_cell_pb.to_csv(OUT_DIR / "pseudobulk_counts_by_celltype.csv", index=False)

    make_pattern_heatmap(global_summary, OUT_DIR / "pattern_gene_directionality_heatmap.png")
    write_readme(OUT_DIR / "README.md", len(genes), global_summary, pattern_summary)
    for suffix, cols in (("global", ["replicate", "condition"]),
                         ("by_celltype", ["replicate", "condition", "cell_type"])):
        group_eligibility(adata.obs, cols, MIN_CELLS_PER_PSEUDOBULK).to_csv(
            OUT_DIR / f"pseudobulk_eligibility_{suffix}.csv", index=False,
        )
    write_method_manifest(OUT_DIR, normalization, CHOSEN_RESULT_H5AD, __file__,
                          MIN_CELLS_PER_PSEUDOBULK, CPM_PSEUDOCOUNT)

    print(f"Wrote outputs to: {OUT_DIR}")
    print(pattern_summary.to_string(index=False))

    if IFN_PATTERN and IFN_PATTERN in set(pattern_summary["pattern"]):
        ifn_pattern = IFN_PATTERN
    else:
        ifn_pattern = (
            pattern_summary
            .sort_values(["n_up_in_stim", "median_mean_log2fc"], ascending=False)
            .iloc[0]["pattern"]
        )
    print(f"\n{ifn_pattern} global gene directions:")
    print(
        global_summary.loc[global_summary["pattern"] == ifn_pattern, [
            "pattern",
            "rank",
            "gene",
            "mean_log2fc_stim_vs_ctrl",
            "n_pairs", "n_positive", "n_negative", "n_zero",
            "agreement_with_mean", "donor_support", "normalization_status",
        ]].to_string(index=False)
    )


if __name__ == "__main__":
    main()
