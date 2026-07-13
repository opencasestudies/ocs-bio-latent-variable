#!/usr/bin/env python3
"""Phase 2 candidate-rank interpretation for the CoGAPS HPC sweep.

This script compares selected full CoGAPS result H5AD files for K=4, K=5,
K=6, and the original K=7 baseline. It intentionally uses only result.obs
(gene loadings), result.var (cell activities), and cached cell metadata.
The large X matrix in each H5AD is not needed.
"""

from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULT_ROOT = ROOT / "results_cogaps_singleprocess_hpc"
RUNS_DIR = RESULT_ROOT / "runs"
CACHE_DIR = RESULT_ROOT / "cache"
OUT_DIR = RESULT_ROOT / "phase2_candidate_interpretation"
FIG_DIR = OUT_DIR / "figures"

PREPROCESSED_H5AD = CACHE_DIR / "preprocessed_cells_hvg3000.h5ad"

STIM_LABEL = "stim"
CTRL_LABEL = "ctrl"

CANDIDATES = [
    {"model": "K4_s2_i2000", "K": 4, "seed": 2, "n_iter": 2000},
    {"model": "K5_s1_i2000", "K": 5, "seed": 1, "n_iter": 2000},
    {"model": "K6_s2_i2000", "K": 6, "seed": 2, "n_iter": 2000},
    {"model": "K7_s2_i2000", "K": 7, "seed": 2, "n_iter": 2000},
]

IFN_MARKERS = {
    "ISG15",
    "ISG20",
    "IFI6",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IRF7",
    "MX1",
    "OAS1",
    "OASL",
    "RSAD2",
    "TNFSF10",
    "BST2",
    "PLSCR1",
    "IFIH1",
    "DDX58",
    "IFITM2",
    "IFITM3",
}


def pattern_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns if str(c).startswith("Pattern")]
    return sorted(cols, key=lambda value: int(str(value).replace("Pattern", "")))


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def safe_corr(a: pd.Series, b: pd.Series) -> float:
    val = a.corr(b)
    if pd.isna(val):
        return float("nan")
    return float(val)


def read_result(candidate: dict[str, object]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    k = candidate["K"]
    seed = candidate["seed"]
    n_iter = candidate["n_iter"]
    tag = f"cogaps_K{k}_seed{seed}_iter{n_iter}"
    result_path = RUNS_DIR / f"{tag}.h5ad"
    metrics_path = RUNS_DIR / f"{tag}.metrics.json"
    if not result_path.exists():
        raise FileNotFoundError(result_path)
    if not metrics_path.exists():
        raise FileNotFoundError(metrics_path)

    result = ad.read_h5ad(result_path, backed="r")
    try:
        pats = pattern_columns(result.obs)
        A = result.obs[pats].copy()
        A.index = A.index.astype(str)
        A.index.name = "gene"
        P = result.var[pats].copy()
        P.index = result.var_names.astype(str)
        P.index.name = "cell_barcode"
    finally:
        result.file.close()

    metrics = json.loads(metrics_path.read_text())
    return A, P, metrics


def load_cell_metadata(cell_index: pd.Index) -> pd.DataFrame:
    adata = ad.read_h5ad(PREPROCESSED_H5AD, backed="r")
    try:
        obs = adata.obs.copy()
    finally:
        adata.file.close()
    obs.index = obs.index.astype(str)
    cols = [c for c in ["condition", "cell_type", "replicate", "label"] if c in obs.columns]
    meta = obs.loc[cell_index.astype(str), cols].copy()
    meta.index.name = "cell_barcode"
    return meta


def top_genes_from_A(A: pd.DataFrame, model: str, n: int = 50) -> pd.DataFrame:
    rows = []
    for pattern in pattern_columns(A):
        top = A[pattern].sort_values(ascending=False).head(n)
        for rank, (gene, weight) in enumerate(top.items(), start=1):
            rows.append(
                {
                    "model": model,
                    "pattern": pattern,
                    "rank": rank,
                    "gene": gene,
                    "weight": float(weight),
                }
            )
    return pd.DataFrame(rows)


def summarize_model(
    candidate: dict[str, object],
    A: pd.DataFrame,
    P: pd.DataFrame,
    meta: pd.DataFrame,
    metrics: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    model = str(candidate["model"])
    pats = pattern_columns(A)
    activity = meta.join(P, how="left")
    stim_indicator = (activity["condition"].astype(str) == STIM_LABEL).astype(int)

    top = top_genes_from_A(A, model=model, n=50)
    top_lookup = {
        pattern: top.query("pattern == @pattern").sort_values("rank")["gene"].tolist()
        for pattern in pats
    }

    pat_corr = P.corr()
    summary_rows = []
    for pattern in pats:
        condition_means = activity.groupby("condition", observed=True)[pattern].mean().to_dict()
        celltype_means = activity.groupby("cell_type", observed=True)[pattern].mean().sort_values(ascending=False)
        replicate_means = activity.groupby("replicate", observed=True)[pattern].mean().sort_values(ascending=False)
        other_corrs = pat_corr.loc[pattern, [p for p in pats if p != pattern]].abs()
        top_genes = top_lookup[pattern]
        ifn_marker_count = len(set(top_genes[:25]) & IFN_MARKERS)
        label = "IFN/ISG candidate" if ifn_marker_count >= 6 or safe_corr(P[pattern], stim_indicator) > 0.6 else "non-IFN candidate"
        summary_rows.append(
            {
                "model": model,
                "K": candidate["K"],
                "seed": candidate["seed"],
                "n_iter": candidate["n_iter"],
                "pattern": pattern,
                "corr_with_stim": safe_corr(P[pattern], stim_indicator),
                "mean_activity_ctrl": float(condition_means.get(CTRL_LABEL, float("nan"))),
                "mean_activity_stim": float(condition_means.get(STIM_LABEL, float("nan"))),
                "stim_minus_ctrl": float(condition_means.get(STIM_LABEL, float("nan")))
                - float(condition_means.get(CTRL_LABEL, float("nan"))),
                "dominant_cell_type": str(celltype_means.index[0]),
                "dominant_cell_type_mean_activity": float(celltype_means.iloc[0]),
                "dominant_replicate": str(replicate_means.index[0]),
                "max_abs_within_model_pattern_corr": float(other_corrs.max()) if len(other_corrs) else float("nan"),
                "ifn_marker_count_top25": ifn_marker_count,
                "candidate_label": label,
                "top_genes_top10": ", ".join(top_genes[:10]),
                "top_genes_top25": ", ".join(top_genes[:25]),
                "metrics_ifn_pattern": metrics.get("ifn_pattern"),
                "metrics_ifn_corr": metrics.get("ifn_corr"),
                "metrics_runtime_min": float(metrics.get("runtime_sec", float("nan"))) / 60,
            }
        )

    condition_rows = []
    for pattern in pats:
        grouped = (
            activity.groupby(["cell_type", "condition"], observed=True)[pattern]
            .agg(["mean", "median", "count"])
            .reset_index()
        )
        grouped.insert(0, "pattern", pattern)
        grouped.insert(0, "model", model)
        grouped = grouped.rename(columns={"mean": "mean_activity", "median": "median_activity", "count": "n_cells"})
        condition_rows.append(grouped)

    replicate_rows = []
    for pattern in pats:
        grouped = (
            activity.groupby(["replicate", "condition"], observed=True)[pattern]
            .agg(["mean", "median", "count"])
            .reset_index()
        )
        grouped.insert(0, "pattern", pattern)
        grouped.insert(0, "model", model)
        grouped = grouped.rename(columns={"mean": "mean_activity", "median": "median_activity", "count": "n_cells"})
        replicate_rows.append(grouped)

    corr_long = (
        pat_corr.reset_index()
        .melt(id_vars="index", var_name="pattern_b", value_name="activity_corr")
        .rename(columns={"index": "pattern_a"})
    )
    corr_long.insert(0, "model", model)

    return (
        pd.DataFrame(summary_rows),
        top,
        pd.concat(condition_rows, ignore_index=True),
        pd.concat(replicate_rows, ignore_index=True),
        corr_long,
    )


def cross_model_tables(
    model_data: dict[str, dict[str, pd.DataFrame]], top_genes: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_sets = {}
    for (model, pattern), group in top_genes.groupby(["model", "pattern"]):
        gene_sets[(model, pattern)] = set(group.sort_values("rank").head(50)["gene"].astype(str))

    top_rows = []
    keys = list(gene_sets)
    for (model_a, pat_a), (model_b, pat_b) in combinations(keys, 2):
        top_rows.append(
            {
                "model_a": model_a,
                "pattern_a": pat_a,
                "model_b": model_b,
                "pattern_b": pat_b,
                "top50_jaccard": jaccard(gene_sets[(model_a, pat_a)], gene_sets[(model_b, pat_b)]),
                "top50_overlap_count": len(gene_sets[(model_a, pat_a)] & gene_sets[(model_b, pat_b)]),
            }
        )
    top_overlap = pd.DataFrame(top_rows)

    activity_rows = []
    model_keys = list(model_data)
    for i, model_a in enumerate(model_keys):
        for model_b in model_keys[i + 1 :]:
            P_a = model_data[model_a]["P"]
            P_b = model_data[model_b]["P"].loc[P_a.index]
            for pat_a in pattern_columns(P_a):
                for pat_b in pattern_columns(P_b):
                    activity_rows.append(
                        {
                            "model_a": model_a,
                            "pattern_a": pat_a,
                            "model_b": model_b,
                            "pattern_b": pat_b,
                            "activity_corr": safe_corr(P_a[pat_a], P_b[pat_b]),
                        }
                    )
    activity_overlap = pd.DataFrame(activity_rows)

    baseline = "K7_s2_i2000"
    k7_activity_rows = []
    k7_gene_rows = []
    for model in model_keys:
        if model == baseline:
            continue
        sub_act = activity_overlap.query(
            "(model_a == @model and model_b == @baseline) or (model_a == @baseline and model_b == @model)"
        ).copy()
        for _, row in sub_act.iterrows():
            if row["model_a"] == baseline:
                row_model = row["model_b"]
                row_pattern = row["pattern_b"]
                k7_pattern = row["pattern_a"]
            else:
                row_model = row["model_a"]
                row_pattern = row["pattern_a"]
                k7_pattern = row["pattern_b"]
            k7_activity_rows.append(
                {
                    "model": row_model,
                    "pattern": row_pattern,
                    "k7_pattern": k7_pattern,
                    "activity_corr": float(row["activity_corr"]),
                }
            )

        sub_gene = top_overlap.query(
            "(model_a == @model and model_b == @baseline) or (model_a == @baseline and model_b == @model)"
        ).copy()
        for _, row in sub_gene.iterrows():
            if row["model_a"] == baseline:
                row_model = row["model_b"]
                row_pattern = row["pattern_b"]
                k7_pattern = row["pattern_a"]
            else:
                row_model = row["model_a"]
                row_pattern = row["pattern_a"]
                k7_pattern = row["pattern_b"]
            k7_gene_rows.append(
                {
                    "model": row_model,
                    "pattern": row_pattern,
                    "k7_pattern": k7_pattern,
                    "top50_jaccard": float(row["top50_jaccard"]),
                    "top50_overlap_count": int(row["top50_overlap_count"]),
                }
            )

    k7_activity = pd.DataFrame(k7_activity_rows)
    k7_gene = pd.DataFrame(k7_gene_rows)

    best_activity = (
        k7_activity.sort_values(["model", "pattern", "activity_corr"], ascending=[True, True, False])
        .groupby(["model", "pattern"], as_index=False)
        .first()
    )
    best_gene = (
        k7_gene.sort_values(["model", "pattern", "top50_jaccard"], ascending=[True, True, False])
        .groupby(["model", "pattern"], as_index=False)
        .first()
    )
    return top_overlap, activity_overlap, best_activity, best_gene


def save_heatmap(
    data: pd.DataFrame,
    outpath: Path,
    title: str,
    cbar_label: str,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    figsize: tuple[float, float] | None = None,
) -> None:
    if figsize is None:
        figsize = (max(7, 0.45 * len(data.columns) + 4), max(4, 0.35 * len(data.index) + 2.5))
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(data.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(data.columns)))
    ax.set_xticklabels(data.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(data.index)))
    ax.set_yticklabels(data.index)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def make_figures(
    pattern_summary: pd.DataFrame,
    celltype_condition: pd.DataFrame,
    best_activity: pd.DataFrame,
    best_gene: pd.DataFrame,
    phase1b_summary: pd.DataFrame,
) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Candidate IFN metric summary.
    ifn_rows = pattern_summary.sort_values("corr_with_stim", ascending=False).groupby("model").first().reset_index()
    order = ["K4_s2_i2000", "K5_s1_i2000", "K6_s2_i2000", "K7_s2_i2000"]
    ifn_rows["model"] = pd.Categorical(ifn_rows["model"], categories=order, ordered=True)
    ifn_rows = ifn_rows.sort_values("model")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(ifn_rows["model"].astype(str), ifn_rows["corr_with_stim"], color=["#4C78A8", "#59A14F", "#F28E2B", "#E15759"])
    ax.set_ylim(0, max(0.8, ifn_rows["corr_with_stim"].max() + 0.03))
    ax.set_ylabel("IFN pattern correlation with stim")
    ax.set_title("Selected candidate IFN pattern strength")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ifn_pattern_correlation_by_candidate.png", dpi=180)
    plt.close(fig)

    # All pattern correlations with stimulation.
    corr_mat = (
        pattern_summary.assign(row=lambda d: d["model"] + ":" + d["pattern"])
        .set_index("row")[["corr_with_stim"]]
        .sort_values("corr_with_stim", ascending=False)
    )
    save_heatmap(
        corr_mat,
        FIG_DIR / "all_pattern_stim_correlations.png",
        "All candidate pattern correlations with stimulation",
        "Pearson r",
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        figsize=(5.5, max(5, 0.28 * len(corr_mat) + 2.5)),
    )

    # Mean pattern activity by condition.
    cond_mat = pattern_summary.assign(row=lambda d: d["model"] + ":" + d["pattern"]).set_index("row")[
        ["mean_activity_ctrl", "mean_activity_stim", "stim_minus_ctrl"]
    ]
    save_heatmap(
        cond_mat,
        FIG_DIR / "pattern_activity_condition_means.png",
        "Pattern activity by condition",
        "activity",
        cmap="viridis",
        figsize=(6, max(5, 0.28 * len(cond_mat) + 2.5)),
    )

    # IFN stim-minus-ctrl by cell type across models.
    ifn_patterns = pattern_summary.sort_values("corr_with_stim", ascending=False).groupby("model").first()["pattern"].to_dict()
    ifn_ct_rows = []
    for model, pattern in ifn_patterns.items():
        sub = celltype_condition.query("model == @model and pattern == @pattern")
        piv = sub.pivot(index="cell_type", columns="condition", values="mean_activity")
        for cell_type, vals in piv.iterrows():
            ifn_ct_rows.append(
                {
                    "model": model,
                    "cell_type": cell_type,
                    "stim_minus_ctrl": vals.get(STIM_LABEL, float("nan")) - vals.get(CTRL_LABEL, float("nan")),
                }
            )
    ifn_ct = pd.DataFrame(ifn_ct_rows)
    ifn_ct_mat = ifn_ct.pivot(index="cell_type", columns="model", values="stim_minus_ctrl")
    ifn_ct_mat = ifn_ct_mat[[c for c in order if c in ifn_ct_mat.columns]]
    save_heatmap(
        ifn_ct_mat,
        FIG_DIR / "ifn_stim_minus_ctrl_by_cell_type.png",
        "IFN pattern stim-minus-control activity by cell type",
        "stim - ctrl",
        cmap="viridis",
        figsize=(7, max(4, 0.4 * len(ifn_ct_mat) + 2)),
    )

    # Best K7 activity matches.
    best_act_mat = best_activity.pivot(index=["model", "pattern"], columns="k7_pattern", values="activity_corr").fillna(0)
    best_act_mat.index = [f"{m}:{p}" for m, p in best_act_mat.index]
    save_heatmap(
        best_act_mat,
        FIG_DIR / "candidate_to_k7_activity_match.png",
        "Candidate pattern activity correlation to K7 patterns",
        "Pearson r",
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        figsize=(8, max(4, 0.35 * len(best_act_mat) + 2.5)),
    )

    # Best K7 top-gene matches.
    best_gene_mat = best_gene.pivot(index=["model", "pattern"], columns="k7_pattern", values="top50_jaccard").fillna(0)
    best_gene_mat.index = [f"{m}:{p}" for m, p in best_gene_mat.index]
    save_heatmap(
        best_gene_mat,
        FIG_DIR / "candidate_to_k7_topgene_jaccard.png",
        "Candidate pattern top-gene Jaccard to K7 patterns",
        "Jaccard",
        cmap="viridis",
        vmin=0,
        vmax=1,
        figsize=(8, max(4, 0.35 * len(best_gene_mat) + 2.5)),
    )

    # Long-iteration sweep metric plot for context.
    fig, ax = plt.subplots(figsize=(7, 4))
    for k, sub in phase1b_summary[phase1b_summary["K"].isin([4, 5, 6, 7])].groupby("K"):
        sub = sub.sort_values("n_iter")
        ax.plot(sub["n_iter"], sub["ifn_corr_mean"], marker="o", label=f"K={k}")
    ax.set_xscale("log")
    ax.set_xlabel("nIterations")
    ax.set_ylabel("Mean IFN correlation across seeds")
    ax.set_title("IFN signal across rank and iteration depth")
    ax.legend(title="Rank")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ifn_correlation_by_k_and_iterations.png", dpi=180)
    plt.close(fig)


def markdown_table(df: pd.DataFrame, cols: list[str], float_cols: list[str] | None = None) -> str:
    float_cols = float_cols or []
    small = df[cols].copy()
    for c in float_cols:
        if c in small:
            small[c] = small[c].map(lambda v: f"{float(v):.4f}" if pd.notna(v) else "")
    headers = list(small.columns)
    rows = []
    for _, row in small.iterrows():
        values = []
        for header in headers:
            value = row[header]
            if pd.isna(value):
                text = ""
            else:
                text = str(value)
            values.append(text.replace("\n", " ").replace("|", "\\|"))
        rows.append(values)
    header_line = "| " + " | ".join(headers) + " |"
    sep_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_line, sep_line, *body])


def write_report(
    pattern_summary: pd.DataFrame,
    best_activity: pd.DataFrame,
    best_gene: pd.DataFrame,
    phase1b_summary: pd.DataFrame,
) -> None:
    report = OUT_DIR / "phase2_candidate_model_interpretation_report.md"
    ifn_rows = pattern_summary.sort_values("corr_with_stim", ascending=False).groupby("model").first().reset_index()
    ifn_rows = ifn_rows.sort_values("K")

    candidate_models = ["K4_s2_i2000", "K5_s1_i2000", "K6_s2_i2000"]
    k7_pattern_counts = []
    for model in candidate_models:
        sub = best_activity.query("model == @model")
        # Count K7 patterns whose best activity match to this candidate model is reasonably strong.
        inv = []
        for k7_pattern, g in sub.groupby("k7_pattern"):
            best = g.sort_values("activity_corr", ascending=False).iloc[0]
            inv.append((k7_pattern, best["pattern"], float(best["activity_corr"])))
        strong = [x for x in inv if abs(x[2]) >= 0.30]
        k7_pattern_counts.append(
            {
                "model": model,
                "num_k7_patterns_with_abs_activity_match_ge_0.30": len(strong),
                "matched_k7_patterns": ", ".join(f"{k7}->{pat} ({corr:.2f})" for k7, pat, corr in strong),
            }
        )
    k7_match_counts = pd.DataFrame(k7_pattern_counts)

    lines = []
    lines.append("# Phase 2 Candidate Model Interpretation\n\n")
    lines.append("## Purpose\n\n")
    lines.append(
        "Phase 2 asks whether the lower-rank candidates that performed well in the sweep "
        "preserve the relevant biology, or whether they collapse biologically meaningful "
        "non-IFN structure that was separated in the original K=7 model.\n\n"
    )
    lines.append("## Candidate Models\n\n")
    lines.append(
        markdown_table(
            ifn_rows,
            [
                "model",
                "K",
                "seed",
                "n_iter",
                "pattern",
                "corr_with_stim",
                "metrics_runtime_min",
                "ifn_marker_count_top25",
                "top_genes_top10",
            ],
            ["corr_with_stim", "metrics_runtime_min"],
        )
    )
    lines.append("\n\n")

    lines.append("## Sweep Context\n\n")
    lines.append(
        markdown_table(
            phase1b_summary[phase1b_summary["K"].isin([4, 5, 6, 7])].sort_values(["K", "n_iter"]),
            ["K", "n_iter", "n_ok", "runtime_min_mean", "ifn_corr_mean", "ifn_corr_std", "jaccard_mean", "jaccard_min"],
            ["runtime_min_mean", "ifn_corr_mean", "ifn_corr_std", "jaccard_mean", "jaccard_min"],
        )
    )
    lines.append("\n\n")

    lines.append("## Pattern Summaries\n\n")
    lines.append(
        markdown_table(
            pattern_summary.sort_values(["K", "pattern"]),
            [
                "model",
                "pattern",
                "corr_with_stim",
                "stim_minus_ctrl",
                "dominant_cell_type",
                "max_abs_within_model_pattern_corr",
                "ifn_marker_count_top25",
                "top_genes_top10",
            ],
            ["corr_with_stim", "stim_minus_ctrl", "max_abs_within_model_pattern_corr"],
        )
    )
    lines.append("\n\n")

    lines.append("## Candidate-to-K7 Matching\n\n")
    lines.append("Best activity match from each candidate pattern to a K=7 pattern:\n\n")
    lines.append(
        markdown_table(
            best_activity.sort_values(["model", "pattern"]),
            ["model", "pattern", "k7_pattern", "activity_corr"],
            ["activity_corr"],
        )
    )
    lines.append("\n\nBest top-gene match from each candidate pattern to a K=7 pattern:\n\n")
    lines.append(
        markdown_table(
            best_gene.sort_values(["model", "pattern"]),
            ["model", "pattern", "k7_pattern", "top50_jaccard", "top50_overlap_count"],
            ["top50_jaccard"],
        )
    )
    lines.append("\n\nK=7 pattern coverage by candidate activity matches:\n\n")
    lines.append(
        markdown_table(
            k7_match_counts,
            [
                "model",
                "num_k7_patterns_with_abs_activity_match_ge_0.30",
                "matched_k7_patterns",
            ],
        )
    )
    lines.append("\n\n")

    lines.append("## Interpretation\n\n")
    lines.append(
        "- `K=4` is the strongest lower-rank candidate by IFN correlation and runtime. "
        "Its IFN-associated pattern is highly coherent and reproducible across the sweep, "
        "but Phase 2 should be read as a check for whether the non-IFN structure is compressed.\n"
    )
    lines.append(
        "- `K=5` and `K=6` recover essentially the same IFN top-gene program while allowing "
        "more non-IFN structure than K=4. Their sweep stability is excellent, especially K=6.\n"
    )
    lines.append(
        "- The K=7 baseline remains useful as a reference model, but the extended sweep no longer "
        "supports K=7 solely on rank-selection metrics. A final decision should favor the lowest K "
        "that preserves interpretable non-IFN programs needed for the narrative.\n"
    )
    lines.append("\n## Phase 2 Recommendation\n\n")
    lines.append(
        "Based on the full pattern-structure comparison, `K=6, seed=2, n_iter=2000` is the strongest "
        "candidate for replacing the original K=7 selected model. `K=4` has the strongest IFN "
        "correlation and best runtime, but it compresses the non-IFN structure: by activity matching "
        "it recovers four K=7 reference programs. `K=5` recovers five K=7 reference programs. `K=6` "
        "recovers six K=7 reference programs, including the main IFN/ISG program, the secondary "
        "monocyte-associated IFN program, monocyte identity/activity structure, CD4/T-cell-associated "
        "structure, and the dendritic/HLA program.\n\n"
    )
    lines.append(
        "The main program not retained as a distinct K=6 factor is K7 Pattern3, a weaker CD8/T-cell-"
        "associated program with mixed housekeeping/cytoskeletal top genes and only mild condition "
        "association. Unless that specific pattern is important to the biological question, "
        "K=6 appears to be the better compromise: lower rank than K=7, near-identical IFN top genes, "
        "excellent seed stability, and better preservation of non-IFN biology than K=4 or K=5.\n\n"
    )
    lines.append(
        "Practical recommendation: use K=6 as the leading selected-rank candidate, keep K=7 as a "
        "reference/sensitivity model, and mention K=4 as the IFN-optimized but biologically compressed "
        "alternative.\n"
    )
    lines.append("\n## Figures\n\n")
    for fig_path in sorted(FIG_DIR.glob("*.png")):
        rel = fig_path.relative_to(OUT_DIR)
        lines.append(f"- `{rel}`\n")

    report.write_text("".join(lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    model_data: dict[str, dict[str, pd.DataFrame]] = {}
    summaries = []
    top_tables = []
    celltype_tables = []
    replicate_tables = []
    within_corr_tables = []

    first_A, first_P, _ = read_result(CANDIDATES[0])
    meta = load_cell_metadata(first_P.index)
    del first_A, first_P

    for candidate in CANDIDATES:
        A, P, metrics = read_result(candidate)
        model = str(candidate["model"])
        model_data[model] = {"A": A, "P": P}
        summary, top, celltype_condition, replicate_condition, within_corr = summarize_model(candidate, A, P, meta, metrics)
        summaries.append(summary)
        top_tables.append(top)
        celltype_tables.append(celltype_condition)
        replicate_tables.append(replicate_condition)
        within_corr_tables.append(within_corr)

    pattern_summary = pd.concat(summaries, ignore_index=True)
    top_genes = pd.concat(top_tables, ignore_index=True)
    celltype_condition = pd.concat(celltype_tables, ignore_index=True)
    replicate_condition = pd.concat(replicate_tables, ignore_index=True)
    within_corr = pd.concat(within_corr_tables, ignore_index=True)
    top_overlap, activity_overlap, best_activity, best_gene = cross_model_tables(model_data, top_genes)

    phase1b_summary = pd.read_csv(RESULT_ROOT / "phase1b_lowerk_long_summary.csv")

    pattern_summary.to_csv(OUT_DIR / "candidate_pattern_summary.csv", index=False)
    top_genes.to_csv(OUT_DIR / "candidate_top_genes.csv", index=False)
    celltype_condition.to_csv(OUT_DIR / "candidate_activity_by_celltype_condition.csv", index=False)
    replicate_condition.to_csv(OUT_DIR / "candidate_activity_by_replicate_condition.csv", index=False)
    within_corr.to_csv(OUT_DIR / "candidate_within_model_activity_correlations.csv", index=False)
    top_overlap.to_csv(OUT_DIR / "candidate_cross_model_topgene_jaccard.csv", index=False)
    activity_overlap.to_csv(OUT_DIR / "candidate_cross_model_activity_correlations.csv", index=False)
    best_activity.to_csv(OUT_DIR / "candidate_best_k7_activity_matches.csv", index=False)
    best_gene.to_csv(OUT_DIR / "candidate_best_k7_topgene_matches.csv", index=False)

    make_figures(pattern_summary, celltype_condition, best_activity, best_gene, phase1b_summary)
    write_report(pattern_summary, best_activity, best_gene, phase1b_summary)

    print(f"Wrote Phase 2 outputs to {OUT_DIR}")
    print(pattern_summary.sort_values("corr_with_stim", ascending=False).groupby("model").first()[["pattern", "corr_with_stim", "top_genes_top10"]])


if __name__ == "__main__":
    main()
