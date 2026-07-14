#!/usr/bin/env python3
"""High-K CoGAPS fragmentation/redundancy diagnostics.

This is a post hoc, non-model-fitting analysis for Phase 3. It reads saved
CoGAPS result H5AD files, uses result.obs as gene loadings (A) and result.var
as cell activities (P), and summarizes whether high-K models produce redundant
or split patterns relative to the lower-rank K6 candidate.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULT_ROOT = ROOT / "results_cogaps_singleprocess_hpc"
RUNS_DIR = RESULT_ROOT / "runs"
OUT_DIR = RESULT_ROOT / "phase3_highk_fragmentation"
JOBS_TSV = RESULT_ROOT / "jobs_phase3_highk_2000.tsv"
K6_REFERENCE = {"K": 6, "seed": 2, "n_iter": 2000}
STIM_LABEL = "stim"
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
    "CXCL10",
    "IFITM1",
    "GBP1",
}


def pattern_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns if str(c).startswith("Pattern")]

    def key(value: str) -> int:
        suffix = str(value).replace("Pattern", "")
        return int(suffix) if suffix.isdigit() else 10**9

    return sorted(cols, key=key)


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def read_result(k: int, seed: int, n_iter: int) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
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


def top_gene_lookup(A: pd.DataFrame, n: int = 50) -> dict[str, list[str]]:
    return {
        pattern: A[pattern].sort_values(ascending=False).head(n).index.astype(str).tolist()
        for pattern in pattern_columns(A)
    }


def offdiag_values(mat: pd.DataFrame, absolute: bool = True) -> list[float]:
    vals = []
    cols = list(mat.columns)
    for a, b in combinations(cols, 2):
        val = mat.loc[a, b]
        if pd.notna(val):
            vals.append(abs(float(val)) if absolute else float(val))
    return vals


def safe_corr(series: pd.Series, target: pd.Series) -> float:
    val = series.corr(target)
    return float(val) if pd.notna(val) else float("nan")


def summarize_one(k: int, seed: int, n_iter: int, ref_ifn_top50: set[str]) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    A, P, metrics = read_result(k, seed, n_iter)
    pats = pattern_columns(A)
    top50 = top_gene_lookup(A, n=50)
    top25 = {p: genes[:25] for p, genes in top50.items()}

    activity_corr = P.corr()
    activity_abs = offdiag_values(activity_corr, absolute=True)

    top_jaccard_rows = []
    top_jaccards = []
    for a, b in combinations(pats, 2):
        jac = jaccard(set(top50[a]), set(top50[b]))
        top_jaccards.append(jac)
        if jac >= 0.25:
            top_jaccard_rows.append(
                {
                    "K": k,
                    "seed": seed,
                    "n_iter": n_iter,
                    "pattern_a": a,
                    "pattern_b": b,
                    "top50_jaccard": jac,
                    "activity_corr": float(activity_corr.loc[a, b]),
                    "pattern_a_top10": ";".join(top50[a][:10]),
                    "pattern_b_top10": ";".join(top50[b][:10]),
                }
            )

    # Condition labels are preserved as cells in result.var for the PyCoGAPS sweep.
    if "condition" in P.index.names:
        raise RuntimeError("Unexpected condition storage")
    condition = None
    # CoGAPS result.var comes from the input cells-by-metadata, so condition is
    # generally available as a column in the full AnnData var dataframe. If it is
    # absent, fall back to marker-count and K6-gene-set based IFN calls.
    result = ad.read_h5ad(RUNS_DIR / f"cogaps_K{k}_seed{seed}_iter{n_iter}.h5ad", backed="r")
    try:
        var_meta = result.var[[c for c in ["condition", "label", "stim"] if c in result.var.columns]].copy()
    finally:
        result.file.close()
    if "condition" in var_meta.columns:
        condition = (var_meta["condition"].astype(str) == STIM_LABEL).astype(int)
    elif "label" in var_meta.columns:
        condition = (var_meta["label"].astype(str) == STIM_LABEL).astype(int)
    elif "stim" in var_meta.columns:
        condition = var_meta["stim"].astype(int)

    pattern_rows = []
    for pattern in pats:
        genes25 = set(top25[pattern])
        genes50 = set(top50[pattern])
        ifn_marker_count_top25 = len(genes25 & IFN_MARKERS)
        jac_to_k6_ifn = jaccard(genes50, ref_ifn_top50)
        corr_with_stim = safe_corr(P[pattern], condition) if condition is not None else float("nan")
        pattern_rows.append(
            {
                "K": k,
                "seed": seed,
                "n_iter": n_iter,
                "pattern": pattern,
                "corr_with_stim": corr_with_stim,
                "ifn_marker_count_top25": ifn_marker_count_top25,
                "top50_jaccard_to_K6_IFN": jac_to_k6_ifn,
                "is_ifn_like": bool(
                    ifn_marker_count_top25 >= 6
                    or jac_to_k6_ifn >= 0.35
                    or (pd.notna(corr_with_stim) and corr_with_stim >= 0.5)
                ),
                "top10_genes": ";".join(top50[pattern][:10]),
            }
        )

    pattern_df = pd.DataFrame(pattern_rows)
    top_pair_df = pd.DataFrame(top_jaccard_rows)
    ifn_like = pattern_df[pattern_df["is_ifn_like"]]
    ifn_metrics_pattern = metrics.get("ifn_pattern")
    ifn_jaccard_to_k6 = float(
        pattern_df.loc[pattern_df["pattern"] == ifn_metrics_pattern, "top50_jaccard_to_K6_IFN"].iloc[0]
    ) if ifn_metrics_pattern in set(pattern_df["pattern"]) else float("nan")

    summary = {
        "K": k,
        "seed": seed,
        "n_iter": n_iter,
        "n_patterns": len(pats),
        "runtime_min": float(metrics.get("runtime_sec", float("nan"))) / 60,
        "metrics_ifn_pattern": ifn_metrics_pattern,
        "metrics_ifn_corr": metrics.get("ifn_corr"),
        "metrics_ifn_top50_jaccard_to_K6_IFN": ifn_jaccard_to_k6,
        "n_ifn_like_patterns": int(len(ifn_like)),
        "max_pattern_jaccard_to_K6_IFN": float(pattern_df["top50_jaccard_to_K6_IFN"].max()),
        "mean_pattern_jaccard_to_K6_IFN": float(pattern_df["top50_jaccard_to_K6_IFN"].mean()),
        "max_abs_activity_corr": max(activity_abs) if activity_abs else float("nan"),
        "mean_abs_activity_corr": float(np.mean(activity_abs)) if activity_abs else float("nan"),
        "n_pattern_pairs_abs_activity_corr_ge_0_8": int(sum(v >= 0.8 for v in activity_abs)),
        "n_pattern_pairs_abs_activity_corr_ge_0_9": int(sum(v >= 0.9 for v in activity_abs)),
        "max_top50_jaccard_between_patterns": max(top_jaccards) if top_jaccards else float("nan"),
        "mean_top50_jaccard_between_patterns": float(np.mean(top_jaccards)) if top_jaccards else float("nan"),
        "n_pattern_pairs_top50_jaccard_ge_0_35": int(sum(v >= 0.35 for v in top_jaccards)),
        "n_pattern_pairs_top50_jaccard_ge_0_50": int(sum(v >= 0.50 for v in top_jaccards)),
    }
    return summary, pattern_df, top_pair_df


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ref_A, _, ref_metrics = read_result(
        int(K6_REFERENCE["K"]),
        int(K6_REFERENCE["seed"]),
        int(K6_REFERENCE["n_iter"]),
    )
    ref_top50 = top_gene_lookup(ref_A, n=50)
    ref_ifn_pattern = ref_metrics.get("ifn_pattern")
    if not ref_ifn_pattern or ref_ifn_pattern not in ref_top50:
        raise RuntimeError(f"Could not determine K6 IFN pattern from metrics: {ref_ifn_pattern}")
    ref_ifn_top50 = set(ref_top50[ref_ifn_pattern])

    jobs = []
    for line in JOBS_TSV.read_text().splitlines():
        if not line.strip():
            continue
        k, seed, n_iter = [int(x) for x in line.split("\t")]
        jobs.append((k, seed, n_iter))

    summary_rows = []
    pattern_tables = []
    pair_tables = []
    for k, seed, n_iter in jobs:
        summary, pattern_df, top_pair_df = summarize_one(k, seed, n_iter, ref_ifn_top50)
        summary_rows.append(summary)
        pattern_tables.append(pattern_df)
        pair_tables.append(top_pair_df)

    summary_df = pd.DataFrame(summary_rows).sort_values(["K", "seed"])
    pattern_df = pd.concat(pattern_tables, ignore_index=True)
    pair_df = pd.concat(pair_tables, ignore_index=True) if pair_tables else pd.DataFrame()

    by_k = (
        summary_df.groupby("K", as_index=False)
        .agg(
            n_ok=("seed", "count"),
            runtime_min_mean=("runtime_min", "mean"),
            metrics_ifn_corr_mean=("metrics_ifn_corr", "mean"),
            metrics_ifn_corr_std=("metrics_ifn_corr", "std"),
            metrics_ifn_top50_jaccard_to_K6_IFN_mean=("metrics_ifn_top50_jaccard_to_K6_IFN", "mean"),
            metrics_ifn_top50_jaccard_to_K6_IFN_min=("metrics_ifn_top50_jaccard_to_K6_IFN", "min"),
            n_ifn_like_patterns_mean=("n_ifn_like_patterns", "mean"),
            n_ifn_like_patterns_max=("n_ifn_like_patterns", "max"),
            max_abs_activity_corr_mean=("max_abs_activity_corr", "mean"),
            max_abs_activity_corr_max=("max_abs_activity_corr", "max"),
            mean_abs_activity_corr_mean=("mean_abs_activity_corr", "mean"),
            max_top50_jaccard_between_patterns_mean=("max_top50_jaccard_between_patterns", "mean"),
            max_top50_jaccard_between_patterns_max=("max_top50_jaccard_between_patterns", "max"),
            n_pattern_pairs_top50_jaccard_ge_0_35_mean=("n_pattern_pairs_top50_jaccard_ge_0_35", "mean"),
            n_pattern_pairs_abs_activity_corr_ge_0_8_mean=("n_pattern_pairs_abs_activity_corr_ge_0_8", "mean"),
        )
        .sort_values("K")
    )

    summary_df.to_csv(OUT_DIR / "phase3_highk_fragmentation_by_run.csv", index=False)
    by_k.to_csv(OUT_DIR / "phase3_highk_fragmentation_by_K.csv", index=False)
    pattern_df.to_csv(OUT_DIR / "phase3_highk_pattern_ifn_fragmentation.csv", index=False)
    pair_df.to_csv(OUT_DIR / "phase3_highk_redundant_pattern_pairs.csv", index=False)
    metadata = {
        "reference_model": K6_REFERENCE,
        "reference_ifn_pattern": ref_ifn_pattern,
        "reference_ifn_top50": sorted(ref_ifn_top50),
        "interpretation": (
            "Counts and pair metrics are descriptive diagnostics of high-K "
            "fragmentation/redundancy. They are not additional model fits."
        ),
    }
    (OUT_DIR / "phase3_highk_fragmentation_metadata.json").write_text(json.dumps(metadata, indent=2))
    print(f"Wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
