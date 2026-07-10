#!/usr/bin/env python3
"""Build learner-facing figures and source tables for selected K6 CoGAPS outputs.

The case-study narrative uses lightweight source tables and saved figure files.
By default, this script creates the R-primary figures and R source tables. It
can also be pointed at the saved Python-compatible output directory to create
parallel Python source tables for the Python tabs.
"""

from __future__ import annotations

import json
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(
    os.environ.get(
        "COGAPS_SELECTED_RESULTS_DIR",
        ROOT / "data/processed/selected_model_k6/r",
    )
)
IMG_DIR = Path(os.environ.get("COGAPS_FIGURES_DIR", ROOT / "data/processed/figures"))
SOURCE_TABLE_DIR = Path(
    os.environ.get(
        "COGAPS_FIGURE_SOURCE_TABLE_DIR",
        ROOT / "data/processed/figures/source_tables/r",
    )
)
LEGACY_IMG_DIR_ENV = os.environ.get("COGAPS_LEGACY_FIGURES_DIR", str(ROOT / "img"))
LEGACY_IMG_DIR = Path(LEGACY_IMG_DIR_ENV) if LEGACY_IMG_DIR_ENV else None
METRICS_JSON = RESULTS_DIR / "cogaps_K6_seed2_iter2000.metrics.json"
PATTERN_SUMMARY = RESULTS_DIR / "pattern_summary.csv"
CELLTYPE_CONDITION = RESULTS_DIR / "pattern_activity_by_celltype_condition.csv"
REPLICATE_CONDITION = RESULTS_DIR / "pattern_activity_by_replicate_condition.csv"
TOP_GENES = RESULTS_DIR / "pattern_top_genes.csv"
DIRECTION_GLOBAL = RESULTS_DIR / "pattern_gene_directionality_global.csv"
CELL_ACTIVITIES = RESULTS_DIR / "pattern_cell_activities_with_metadata.csv.gz"
DIRECTIONALITY_HEATMAP = RESULTS_DIR / "pattern_gene_directionality_heatmap.png"
PSEUDOBULK_GLOBAL = RESULTS_DIR / "pseudobulk_counts_global.csv"
PSEUDOBULK_BY_CELLTYPE = RESULTS_DIR / "pseudobulk_counts_by_celltype.csv"

CONTROL_COLOR = "#64748B"
STIM_COLOR = "#E76F51"
IFN_COLOR = "#E76F51"
OTHER_COLOR = "#64748B"
BLUE = "#2B8CBE"
TEAL = "#168477"
GRID = "#E2E8F0"
TEXT = "#0B2239"
CELL_TYPE_COLORS = {
    "CD4 T cells": "#2A9D8F",
    "CD14+ Monocytes": "#457B9D",
    "B cells": "#E9C46A",
    "NK cells": "#7E68B1",
    "CD8 T cells": "#52B788",
    "FCGR3A+ Monocytes": "#F4A261",
    "Dendritic cells": "#5B6770",
    "Megakaryocytes": "#D1495B",
}

PATTERN_WORKING_LABELS = {
    "Pattern1": (
        "FCGR3A+ monocyte-linked",
        "Activity is strongest in the FCGR3A+ monocyte neighborhood; stimulation association is weak.",
    ),
    "Pattern2": (
        "CD14+ monocyte-linked",
        "Activity is strongest in CD14+ monocytes and is lower after IFN-beta stimulation.",
    ),
    "Pattern3": (
        "Dendritic / antigen-presentation-linked",
        "High-weight HLA genes and dendritic-cell activity point to antigen-presentation structure.",
    ),
    "Pattern4": (
        "IFN-associated stimulation",
        "Canonical interferon-stimulated genes and strong positive stimulation correlation mark the main IFN response.",
    ),
    "Pattern5": (
        "candidate secondary IFN / myeloid response",
        "Positive stimulation association and CXCL10/APOBEC3A genes suggest a second myeloid response signal.",
    ),
    "Pattern6": (
        "T-cell / baseline lymphocyte-linked",
        "Activity is strongest in T-cell neighborhoods and is lower after stimulation.",
    ),
}

plt.rcParams.update(
    {
        "font.size": 12,
        "axes.titlesize": 18,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
        "axes.edgecolor": "#CBD5E1",
        "axes.labelcolor": TEXT,
        "xtick.color": "#334155",
        "ytick.color": "#334155",
        "text.color": TEXT,
    }
)


def pattern_working_interpretation(
    pattern: str,
    dominant_cell_type: str,
    top_genes: str = "",
    corr_with_stim: float | None = None,
    ifn_pattern: str | None = None,
) -> tuple[str, str]:
    """Return cautious learner-facing labels for selected-model patterns."""
    if ifn_pattern is not None and str(pattern) == str(ifn_pattern):
        return (
            "IFN-associated stimulation",
            "Canonical interferon-stimulated genes and strong positive stimulation correlation mark the main IFN response.",
        )

    genes = str(top_genes).upper()
    if corr_with_stim is not None and corr_with_stim > 0.25:
        if any(marker in genes for marker in ["ISG", "IFIT", "MX1", "OAS", "IRF7", "RSAD2", "CXCL10", "APOBEC3A"]):
            return (
                "candidate secondary IFN / myeloid response",
                "Positive stimulation association and interferon-related genes suggest a second response-linked signal.",
            )

    if "Dendritic" in str(dominant_cell_type) or "HLA-" in genes or "CD74" in genes:
        return (
            "Dendritic / antigen-presentation-linked",
            "High-weight HLA genes and dendritic-cell activity point to antigen-presentation structure.",
        )
    if "FCGR3A" in str(dominant_cell_type) or "FCGR3A" in genes:
        return (
            "FCGR3A+ monocyte-linked",
            "Activity is strongest in the FCGR3A+ monocyte neighborhood; stimulation association is weak.",
        )
    if "CD14" in str(dominant_cell_type):
        return (
            "CD14+ monocyte-linked",
            "Activity is strongest in CD14+ monocytes and is lower after IFN-beta stimulation.",
        )
    if any(label in str(dominant_cell_type) for label in ["CD4 T cells", "CD8 T cells", "NK cells"]):
        return (
            "T-cell / baseline lymphocyte-linked",
            "Activity is strongest in lymphocyte neighborhoods and is lower after stimulation.",
        )
    if "B cells" in str(dominant_cell_type):
        return (
            "B-cell-linked",
            "Activity is strongest in the B-cell neighborhood.",
        )
    return (
        f"{dominant_cell_type}-linked",
        f"Activity is most concentrated in the {dominant_cell_type} neighborhood.",
    )


def pattern_columns(df: pd.DataFrame) -> list[str]:
    return sorted(
        [c for c in df.columns if str(c).startswith("Pattern")],
        key=lambda value: int(str(value).replace("Pattern", "")),
    )


def pattern_sort_key(value: str) -> int:
    return int(str(value).replace("Pattern", ""))


def sorted_patterns(patterns: list[str] | pd.Series) -> list[str]:
    return sorted([str(p) for p in patterns], key=pattern_sort_key)


def selected_ifn_pattern() -> str:
    with METRICS_JSON.open() as f:
        metrics = json.load(f)
    return str(metrics["ifn_pattern"])


def save_figure(fig: plt.Figure, filename: str, dpi: int = 220) -> None:
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    if LEGACY_IMG_DIR is not None:
        LEGACY_IMG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMG_DIR / filename, dpi=dpi, bbox_inches="tight")
    if LEGACY_IMG_DIR is not None:
        fig.savefig(LEGACY_IMG_DIR / filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def save_source_table(df: pd.DataFrame, filename: str) -> None:
    SOURCE_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(SOURCE_TABLE_DIR / filename, index=False)


def save_source_table_gz(df: pd.DataFrame, filename: str) -> None:
    SOURCE_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(SOURCE_TABLE_DIR / filename, index=False, compression="gzip")


def style_axis(ax: plt.Axes, xgrid: bool = False, ygrid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x" if xgrid else "y", color=GRID, linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        color="white",
        ha="center",
        va="center",
        bbox=dict(boxstyle="round,pad=0.28", facecolor=TEAL, edgecolor="none"),
    )


def save_pattern_condition_overview(ifn_pattern: str) -> None:
    summary = pd.read_csv(PATTERN_SUMMARY)
    summary["stim_minus_ctrl"] = summary["mean_activity_stim"] - summary["mean_activity_ctrl"]
    summary["is_ifn_candidate"] = summary["pattern"].eq(ifn_pattern)
    interpretations = summary.apply(
        lambda row: pattern_working_interpretation(
            row["pattern"],
            row["dominant_cell_type"],
            row.get("top_genes_top8", ""),
            row.get("corr_with_stim"),
            ifn_pattern,
        ),
        axis=1,
    )
    summary["working_label"] = [item[0] for item in interpretations]
    summary["interpretation_note"] = [item[1] for item in interpretations]
    summary = summary.sort_values("pattern", key=lambda s: s.map(pattern_sort_key))

    source = summary[
        [
            "pattern",
            "working_label",
            "interpretation_note",
            "corr_with_stim",
            "mean_activity_ctrl",
            "mean_activity_stim",
            "stim_minus_ctrl",
            "dominant_cell_type",
            "top_genes_top8",
            "is_ifn_candidate",
        ]
    ].copy()
    save_source_table(source, "figure_viz_pattern_condition_overview_source.csv")

    patterns = source["pattern"].tolist()
    y = np.arange(len(patterns))
    colors = np.where(source["is_ifn_candidate"], IFN_COLOR, OTHER_COLOR)
    y_labels = [f"{pattern}\n{label}" for pattern, label in zip(source["pattern"], source["working_label"])]

    fig, axes = plt.subplots(
        ncols=2,
        figsize=(14.6, 6.2),
        gridspec_kw={"width_ratios": [1.15, 1.05], "wspace": 0.30},
    )

    ax = axes[0]
    ax.axvline(0, color="#CBD5E1", linewidth=1.4)
    for yi, corr, color in zip(y, source["corr_with_stim"], colors):
        ax.plot([0, corr], [yi, yi], color=color, linewidth=4, solid_capstyle="round")
    ax.scatter(source["corr_with_stim"], y, s=95, c=colors, edgecolor="white", linewidth=1.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Correlation with IFN-beta stimulation")
    ax.set_title("Which patterns track stimulation?")
    ax.set_xlim(-0.35, 0.82)
    style_axis(ax, xgrid=True, ygrid=False)
    add_panel_label(ax, "A")

    ax = axes[1]
    ax.hlines(y, source["mean_activity_ctrl"], source["mean_activity_stim"], color="#CBD5E1", linewidth=3)
    ax.scatter(source["mean_activity_ctrl"], y, s=82, color=CONTROL_COLOR, edgecolor="white", linewidth=1.2, label="control", zorder=3)
    ax.scatter(source["mean_activity_stim"], y, s=82, color=STIM_COLOR, edgecolor="white", linewidth=1.2, label="IFN-beta, 6 h", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.set_xlabel("Mean pattern activity")
    ax.set_title("How much does activity shift?")
    ax.legend(frameon=False, loc="lower right")
    style_axis(ax, xgrid=True, ygrid=False)
    add_panel_label(ax, "B")

    fig.suptitle("Pattern-level screen for condition-associated activity", fontsize=21, fontweight="bold", y=1.05)
    fig.text(
        0.5,
        -0.02,
        "The highlighted pattern has the strongest positive association with IFN-beta stimulation in the selected R run.",
        ha="center",
        fontsize=12,
        color="#475569",
    )
    save_figure(fig, "figure_viz_pattern_condition_overview.png")


def celltype_centroids(source: pd.DataFrame) -> pd.DataFrame:
    """Return robust label positions for cell-type labels on the UMAP."""
    return (
        source.groupby("cell_type", observed=True)[["umap_1", "umap_2"]]
        .median()
        .reset_index()
    )


def short_celltype_label(label: str) -> str:
    replacements = {
        "CD14+ Monocytes": "CD14+\nmonocytes",
        "FCGR3A+ Monocytes": "FCGR3A+\nmonocytes",
        "Dendritic cells": "Dendritic\ncells",
        "Megakaryocytes": "Megakaryocytes",
        "CD4 T cells": "CD4 T cells",
        "CD8 T cells": "CD8 T cells",
        "B cells": "B cells",
        "NK cells": "NK cells",
    }
    return replacements.get(label, label)


def draw_identity_reference(
    ax: plt.Axes,
    source: pd.DataFrame,
    centroids: pd.DataFrame,
    point_size: float = 2.0,
) -> None:
    for cell_type, group in source.groupby("cell_type", sort=False):
        ax.scatter(
            group["umap_1"],
            group["umap_2"],
            s=point_size,
            alpha=0.75,
            color=CELL_TYPE_COLORS.get(cell_type, "#94A3B8"),
            linewidths=0,
        )
    for _, row in centroids.iterrows():
        ax.text(
            row["umap_1"],
            row["umap_2"],
            short_celltype_label(row["cell_type"]),
            ha="center",
            va="center",
            fontsize=10.5,
            fontweight="bold",
            color=TEXT,
            bbox=dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                edgecolor="#CBD5E1",
                alpha=0.94,
            ),
        )
    ax.set_title("Cell identity reference", fontweight="bold")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")


def save_umap_pattern_pair_figures(df: pd.DataFrame, ifn_pattern: str) -> None:
    """Save one large identity-vs-activity UMAP figure per selected-model pattern."""
    summary = pd.read_csv(PATTERN_SUMMARY)
    patterns = sorted_patterns(summary["pattern"])
    base_cols = ["cell_barcode", "cell_type", "condition", "umap_1", "umap_2"]
    source_wide = df[base_cols + patterns].copy()

    long = source_wide.melt(
        id_vars=base_cols,
        value_vars=patterns,
        var_name="pattern",
        value_name="activity",
    )
    long = long.merge(
        summary[["pattern", "corr_with_stim", "dominant_cell_type", "top_genes_top8"]],
        on="pattern",
        how="left",
    )
    interpretations = long.apply(
        lambda row: pattern_working_interpretation(
            row["pattern"],
            row["dominant_cell_type"],
            row.get("top_genes_top8", ""),
            row.get("corr_with_stim"),
            ifn_pattern,
        ),
        axis=1,
    )
    long["working_label"] = [item[0] for item in interpretations]
    long["interpretation_note"] = [item[1] for item in interpretations]
    long["is_ifn_candidate"] = long["pattern"].eq(ifn_pattern)
    save_source_table_gz(long, "figure_viz_umap_pattern_activity_pairs_source.csv.gz")

    centroids = celltype_centroids(source_wide)
    global_vmax = float(np.nanquantile(source_wide[patterns].to_numpy().ravel(), 0.995))
    global_vmax = max(global_vmax, 0.1)
    x_min, x_max = source_wide["umap_1"].min(), source_wide["umap_1"].max()
    y_min, y_max = source_wide["umap_2"].min(), source_wide["umap_2"].max()
    x_pad = (x_max - x_min) * 0.04
    y_pad = (y_max - y_min) * 0.04

    for pattern in patterns:
        pattern_df = source_wide.sort_values(pattern, ascending=True)
        row = summary.loc[summary["pattern"].eq(pattern)].iloc[0]
        filename = f"figure_viz_umap_{pattern.lower()}_identity_activity.png"
        display_pattern = pattern.replace("Pattern", "Pattern ")
        corr = float(row["corr_with_stim"])
        dominant = str(row["dominant_cell_type"])
        working_label, _ = pattern_working_interpretation(
            pattern,
            dominant,
            row.get("top_genes_top8", ""),
            row.get("corr_with_stim"),
            ifn_pattern,
        )

        title = f"{display_pattern}: {working_label}"

        fig = plt.figure(figsize=(15.2, 7.8))
        grid = fig.add_gridspec(
            nrows=1,
            ncols=3,
            width_ratios=[1.0, 1.0, 0.045],
            wspace=0.07,
        )
        ax_identity = fig.add_subplot(grid[0, 0])
        ax_activity = fig.add_subplot(grid[0, 1])
        cax = fig.add_subplot(grid[0, 2])

        draw_identity_reference(ax_identity, source_wide, centroids, point_size=4.2)

        sc = ax_activity.scatter(
            pattern_df["umap_1"],
            pattern_df["umap_2"],
            c=pattern_df[pattern],
            cmap="magma",
            vmin=0,
            vmax=global_vmax,
            s=4.2,
            alpha=0.86,
            linewidths=0,
        )
        ax_activity.set_title(
            f"{display_pattern} activity",
            fontweight="bold",
        )
        ax_activity.set_xlabel("UMAP 1")
        ax_activity.set_ylabel("UMAP 2")
        cbar = fig.colorbar(sc, cax=cax)
        cbar.set_label("Pattern activity")

        for ax in (ax_identity, ax_activity):
            ax.set_xlim(x_min - x_pad, x_max + x_pad)
            ax.set_ylim(y_min - y_pad, y_max + y_pad)
            ax.set_aspect("equal", adjustable="box")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(False)

        fig.suptitle(title, fontsize=24, fontweight="bold", y=0.98)
        fig.subplots_adjust(top=0.84, bottom=0.12, left=0.055, right=0.95)
        save_figure(fig, filename, dpi=230)


def save_ifn_replicate_pairs(ifn_pattern: str) -> None:
    replicate = pd.read_csv(REPLICATE_CONDITION)
    replicate = replicate.loc[replicate["pattern"].eq(ifn_pattern)].copy()
    wide = (
        replicate.pivot_table(
            index="replicate",
            columns="condition",
            values="mean_activity",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    wide["stim_minus_ctrl"] = wide["stim"] - wide["ctrl"]
    wide = wide.sort_values("stim_minus_ctrl", ascending=False)
    wide["donor_label"] = [f"Donor {i}" for i in range(1, len(wide) + 1)]
    save_source_table(wide, "figure_viz_ifn_replicate_pairs_source.csv")

    positive_n = int((wide["stim_minus_ctrl"] > 0).sum())
    mean_delta = float(wide["stim_minus_ctrl"].mean())

    fig, ax = plt.subplots(figsize=(9.4, 5.9))
    y = np.arange(len(wide))
    ax.axvline(0, color="#CBD5E1", linewidth=1.5)
    ax.hlines(y, 0, wide["stim_minus_ctrl"], color="#F2B7A4", linewidth=7, alpha=0.95)
    ax.scatter(
        wide["stim_minus_ctrl"],
        y,
        color=STIM_COLOR,
        s=92,
        edgecolor="white",
        linewidth=1.3,
        zorder=3,
    )
    ax.axvline(mean_delta, color=TEXT, linewidth=2.8, linestyle="--")
    ax.text(
        mean_delta + 0.004,
        len(wide) - 1.55,
        f"mean shift = {mean_delta:.2f}",
        ha="left",
        va="center",
        fontsize=11,
        color=TEXT,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#CBD5E1"),
    )
    ax.set_yticks(y)
    ax.set_yticklabels(wide["donor_label"])
    ax.invert_yaxis()
    ax.set_xlabel("Mean activity difference (IFN-beta - control)")
    ax.set_ylabel("Matched donor pair")
    ax.set_title("Donor-paired shifts in IFN-associated activity")
    ax.set_xlim(0, max(0.34, wide["stim_minus_ctrl"].max() * 1.15))
    style_axis(ax, xgrid=True, ygrid=False)
    fig.subplots_adjust(top=0.78, bottom=0.16, left=0.18, right=0.96)
    fig.text(
        0.68,
        0.9,
        f"{positive_n} of {len(wide)} donor pairs are higher after stimulation",
        ha="center",
        va="center",
        fontsize=12,
        color=TEXT,
        bbox=dict(boxstyle="round,pad=0.42", facecolor="#EEF6F4", edgecolor="#B8DDD6"),
    )
    fig.text(
        0.5,
        0.035,
        "Each row compares matched control and stimulated aliquots from the same donor.",
        ha="center",
        fontsize=12,
        color="#475569",
    )
    save_figure(fig, "figure_viz_ifn_replicate_pairs.png")


def save_ifn_celltype_condition(ifn_pattern: str) -> None:
    celltype = pd.read_csv(CELLTYPE_CONDITION)
    ifn = celltype.loc[celltype["pattern"].eq(ifn_pattern)].copy()
    wide = (
        ifn.pivot_table(
            index="cell_type",
            columns="condition",
            values="mean_activity",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    counts = (
        ifn.pivot_table(index="cell_type", columns="condition", values="n_cells", aggfunc="first")
        .reset_index()
        .rename_axis(None, axis=1)
        .rename(columns={"ctrl": "n_ctrl", "stim": "n_stim"})
    )
    wide = wide.merge(counts, on="cell_type", how="left")
    wide["stim_minus_ctrl"] = wide["stim"] - wide["ctrl"]
    wide = wide.sort_values("stim", ascending=True)
    save_source_table(wide, "figure_viz_ifn_celltype_condition_source.csv")

    y = np.arange(len(wide))
    fig, ax = plt.subplots(figsize=(9.2, 6.3))
    ax.hlines(y, wide["ctrl"], wide["stim"], color="#CBD5E1", linewidth=3)
    ax.scatter(wide["ctrl"], y, s=92, color=CONTROL_COLOR, edgecolor="white", linewidth=1.2, label="control", zorder=3)
    ax.scatter(wide["stim"], y, s=92, color=STIM_COLOR, edgecolor="white", linewidth=1.2, label="IFN-beta, 6 h", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(wide["cell_type"])
    ax.set_xlabel("Mean IFN-associated pattern activity")
    ax.set_title("IFN-associated activity is broad but not uniform")
    ax.legend(frameon=False, loc="lower right")
    style_axis(ax, xgrid=True, ygrid=False)
    fig.text(
        0.5,
        -0.03,
        "A broad response should rise in several cell types, while different heights remind us that cell identity still matters.",
        ha="center",
        fontsize=12,
        color="#475569",
    )
    save_figure(fig, "figure_viz_ifn_celltype_condition.png")


def save_pattern_celltype_shift_heatmap() -> None:
    celltype = pd.read_csv(CELLTYPE_CONDITION)
    wide = (
        celltype.pivot_table(
            index=["pattern", "cell_type"],
            columns="condition",
            values="mean_activity",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    wide["stim_minus_ctrl"] = wide["stim"] - wide["ctrl"]
    wide = wide.sort_values(["pattern", "cell_type"], key=lambda s: s.map(pattern_sort_key) if s.name == "pattern" else s)
    save_source_table(wide, "figure_viz_pattern_celltype_shift_source.csv")

    patterns = sorted_patterns(wide["pattern"].unique())
    cell_types = (
        wide.groupby("cell_type")["stim_minus_ctrl"]
        .apply(lambda x: float(np.nanmax(np.abs(x))))
        .sort_values(ascending=False)
        .index.tolist()
    )
    matrix = (
        wide.pivot(index="pattern", columns="cell_type", values="stim_minus_ctrl")
        .reindex(index=patterns, columns=cell_types)
    )
    values = matrix.to_numpy()
    vmax = np.nanmax(np.abs(values))
    vmax = max(vmax, 0.01)

    fig, ax = plt.subplots(figsize=(12.2, 6.7))
    im = ax.imshow(values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(patterns)))
    ax.set_yticklabels(patterns)
    ax.set_title("Condition shifts by pattern and PBMC identity")
    ax.set_xlabel("Annotated PBMC identity")
    ax.set_ylabel("CoGAPS pattern")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="#111827")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Mean activity: IFN-beta - control")
    fig.subplots_adjust(left=0.12, right=0.84, bottom=0.32, top=0.88)
    fig.text(
        0.5,
        0.04,
        "Positive values show higher mean activity in stimulated cells; negative values show higher mean activity in controls.",
        ha="center",
        fontsize=12,
        color="#475569",
    )
    save_figure(fig, "figure_viz_pattern_celltype_shift_heatmap.png")


def save_top_gene_directionality(ifn_pattern: str) -> None:
    direction = pd.read_csv(DIRECTION_GLOBAL)
    selected = direction.loc[
        direction["pattern"].eq(ifn_pattern) & direction["rank"].le(15)
    ].copy()
    selected = selected.sort_values("rank", ascending=False)
    save_source_table(selected, "figure_viz_top_gene_directionality_source.csv")

    y = np.arange(len(selected))
    colors = np.where(selected["direction"].eq("up_in_stim"), STIM_COLOR, BLUE)

    fig, axes = plt.subplots(
        ncols=2,
        figsize=(11.5, 6.8),
        sharey=True,
        gridspec_kw={"width_ratios": [0.95, 1.25], "wspace": 0.08},
    )

    ax = axes[0]
    ax.barh(y, selected["cogaps_weight"], color="#94A3B8", edgecolor="white")
    ax.set_yticks(y)
    ax.set_yticklabels(selected["gene"])
    ax.set_xlabel("CoGAPS gene weight")
    ax.set_title("Genes that define the pattern")
    style_axis(ax, xgrid=True, ygrid=False)

    ax = axes[1]
    ax.axvline(0, color="#CBD5E1", linewidth=1.5)
    ax.barh(y, selected["mean_log2fc_stim_vs_ctrl"], color=colors, edgecolor="white")
    ax.set_xlabel("Mean log2 fold-change\nIFN-beta vs control")
    ax.set_title("Direction in stimulation")
    style_axis(ax, xgrid=True, ygrid=False)
    handles = [
        plt.Line2D([0], [0], color=STIM_COLOR, lw=8, label="higher in stimulation"),
        plt.Line2D([0], [0], color=BLUE, lw=8, label="lower in stimulation"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")

    fig.suptitle("Top genes need both weight and direction", fontsize=21, fontweight="bold", y=1.03)
    fig.text(
        0.5,
        -0.02,
        "A high weight says a gene helps define the pattern; the expression contrast says whether it rises or falls after IFN-beta.",
        ha="center",
        fontsize=12,
        color="#475569",
    )
    save_figure(fig, "figure_viz_top_gene_directionality.png")


def selected_identity_pattern(ifn_pattern: str) -> str:
    """Choose a non-IFN pattern that gives a clear identity-linked contrast."""
    summary = pd.read_csv(PATTERN_SUMMARY)
    summary = summary.loc[~summary["pattern"].eq(ifn_pattern)].copy()
    summary["abs_corr_with_stim"] = summary["corr_with_stim"].abs()
    marker_sets = summary["top_genes_top15"].fillna("").str.upper()
    antigen_mask = marker_sets.str.contains("CD74|HLA-DRA|HLA-DRB1|HLA-DPA1|HLA-DPB1", regex=True)
    if antigen_mask.any():
        candidates = summary.loc[antigen_mask].sort_values(
            ["abs_corr_with_stim", "pattern"],
            key=lambda col: col.map(pattern_sort_key) if col.name == "pattern" else col,
        )
        return str(candidates.iloc[0]["pattern"])

    candidates = summary.sort_values(
        ["abs_corr_with_stim", "pattern"],
        key=lambda col: col.map(pattern_sort_key) if col.name == "pattern" else col,
    )
    return str(candidates.iloc[0]["pattern"])


def cpm_column_for_gene(gene: str, df: pd.DataFrame) -> str | None:
    """Find the CPM column for a gene, handling R/Python name punctuation."""
    gene = str(gene)
    candidates = [
        f"{gene}__cpm",
        f"{gene.replace('.', '-')}__cpm",
        f"{gene.replace('-', '.')}__cpm",
    ]
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    return None


def expression_per_million(gene: str, df: pd.DataFrame) -> pd.Series:
    """Return CPM-like expression for a gene, falling back to raw counts."""
    cpm_col = cpm_column_for_gene(gene, df)
    if cpm_col is not None and df[cpm_col].notna().any():
        return pd.to_numeric(df[cpm_col], errors="coerce")

    gene = str(gene)
    raw_candidates = [
        gene,
        gene.replace(".", "-"),
        gene.replace("-", "."),
    ]
    raw_col = next((candidate for candidate in raw_candidates if candidate in df.columns), None)
    if raw_col is None:
        return pd.Series(np.nan, index=df.index)

    raw = pd.to_numeric(df[raw_col], errors="coerce")
    if "library_size_total_counts" not in df.columns:
        return raw
    library_size = pd.to_numeric(df["library_size_total_counts"], errors="coerce")
    return raw / library_size.replace(0, np.nan) * 1_000_000


def mean_stim_log2fc_from_global_pseudobulk(gene: str, pseudobulk: pd.DataFrame) -> float:
    pseudobulk = pseudobulk.copy()
    pseudobulk["_expr_per_million"] = expression_per_million(gene, pseudobulk)

    values = []
    for _, group in pseudobulk.groupby("replicate", observed=True):
        ctrl = group.loc[group["condition"].eq("ctrl"), "_expr_per_million"].dropna()
        stim = group.loc[group["condition"].eq("stim"), "_expr_per_million"].dropna()
        if ctrl.empty or stim.empty:
            continue
        values.append(np.log2((float(stim.iloc[0]) + 1) / (float(ctrl.iloc[0]) + 1)))
    return float(np.nanmean(values)) if values else np.nan


def mean_identity_log2fc_from_celltype_pseudobulk(
    gene: str,
    dominant_cell_type: str,
    pseudobulk: pd.DataFrame,
) -> float:
    pseudobulk = pseudobulk.copy()
    pseudobulk["_expr_per_million"] = expression_per_million(gene, pseudobulk)

    values = []
    for _, group in pseudobulk.groupby(["replicate", "condition"], observed=True):
        dominant = group.loc[group["cell_type"].eq(dominant_cell_type), "_expr_per_million"].dropna()
        others = group.loc[
            ~group["cell_type"].eq(dominant_cell_type),
            ["_expr_per_million", "n_cells"],
        ].dropna()
        if dominant.empty or others.empty:
            continue
        if "n_cells" in others.columns and others["n_cells"].sum() > 0:
            other_mean = np.average(others["_expr_per_million"], weights=others["n_cells"])
        else:
            other_mean = others["_expr_per_million"].mean()
        values.append(np.log2((float(dominant.iloc[0]) + 1) / (float(other_mean) + 1)))
    return float(np.nanmean(values)) if values else np.nan


def save_identity_pattern_directionality(identity_pattern: str) -> None:
    summary = pd.read_csv(PATTERN_SUMMARY)
    direction = pd.read_csv(DIRECTION_GLOBAL)
    pseudobulk_global = pd.read_csv(PSEUDOBULK_GLOBAL)
    pseudobulk_celltype = pd.read_csv(PSEUDOBULK_BY_CELLTYPE)

    summary_row = summary.loc[summary["pattern"].eq(identity_pattern)].iloc[0]
    dominant_cell_type = str(summary_row["dominant_cell_type"])
    label, note = pattern_working_interpretation(
        identity_pattern,
        dominant_cell_type,
        summary_row.get("top_genes_top8", ""),
        summary_row.get("corr_with_stim"),
        ifn_pattern=None,
    )

    selected = direction.loc[
        direction["pattern"].eq(identity_pattern) & direction["rank"].le(12)
    ].copy()
    selected["dominant_cell_type"] = dominant_cell_type
    selected["working_label"] = label
    selected["interpretation_note"] = note
    selected["corr_with_stim"] = float(summary_row["corr_with_stim"])
    selected["identity_log2fc_dominant_vs_other"] = [
        mean_identity_log2fc_from_celltype_pseudobulk(gene, dominant_cell_type, pseudobulk_celltype)
        for gene in selected["gene"]
    ]
    selected["stim_log2fc_stim_vs_ctrl"] = [
        mean_stim_log2fc_from_global_pseudobulk(gene, pseudobulk_global)
        for gene in selected["gene"]
    ]
    selected["identity_direction"] = np.where(
        selected["identity_log2fc_dominant_vs_other"].ge(0),
        f"higher in {dominant_cell_type}",
        "higher in other PBMCs",
    )
    selected["condition_direction"] = np.where(
        selected["stim_log2fc_stim_vs_ctrl"].ge(0),
        "higher in stimulation",
        "higher in control",
    )
    selected = selected.sort_values("rank", ascending=False)
    save_source_table(selected, "figure_viz_identity_pattern_directionality_source.csv")

    y = np.arange(len(selected))
    identity_colors = np.where(selected["identity_log2fc_dominant_vs_other"].ge(0), TEAL, "#7E68B1")
    condition_colors = np.where(selected["stim_log2fc_stim_vs_ctrl"].ge(0), STIM_COLOR, BLUE)

    fig, axes = plt.subplots(
        ncols=3,
        figsize=(14, 6.7),
        sharey=True,
        gridspec_kw={"width_ratios": [0.9, 1.2, 1.15], "wspace": 0.08},
    )

    axes[0].barh(y, selected["cogaps_weight"], color="#94A3B8", edgecolor="white")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(selected["gene"])
    axes[0].set_xlabel("CoGAPS gene weight")
    axes[0].set_title("Genes that define\nthis pattern")
    style_axis(axes[0], xgrid=True, ygrid=False)

    axes[1].axvline(0, color="#CBD5E1", linewidth=1.5)
    axes[1].barh(
        y,
        selected["identity_log2fc_dominant_vs_other"],
        color=identity_colors,
        edgecolor="white",
    )
    axes[1].set_xlabel(f"Mean log2 fold-change\n{dominant_cell_type} vs other PBMCs")
    axes[1].set_title("Direction across\ncell identities")
    style_axis(axes[1], xgrid=True, ygrid=False)

    axes[2].axvline(0, color="#CBD5E1", linewidth=1.5)
    axes[2].barh(
        y,
        selected["stim_log2fc_stim_vs_ctrl"],
        color=condition_colors,
        edgecolor="white",
    )
    axes[2].set_xlabel("Mean log2 fold-change\nIFN-beta vs control")
    axes[2].set_title("Direction across\nconditions")
    style_axis(axes[2], xgrid=True, ygrid=False)

    handles = [
        plt.Line2D([0], [0], color=TEAL, lw=8, label=f"higher in {dominant_cell_type}"),
        plt.Line2D([0], [0], color="#7E68B1", lw=8, label="higher in other PBMCs"),
        plt.Line2D([0], [0], color=STIM_COLOR, lw=8, label="higher in stimulation"),
        plt.Line2D([0], [0], color=BLUE, lw=8, label="higher in control"),
    ]
    fig.legend(handles=handles, frameon=False, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.045))

    fig.suptitle("Identity-linked top genes need a different direction check", fontsize=20, fontweight="bold", y=1.02)
    fig.text(
        0.5,
        -0.095,
        f"{identity_pattern} is {label}; compare whether its top genes mark {dominant_cell_type}, stimulation, or both.",
        ha="center",
        fontsize=12,
        color="#475569",
    )
    save_figure(fig, "figure_viz_identity_pattern_directionality.png")


def save_ifn_by_celltype_split(df: pd.DataFrame, ifn_pattern: str) -> None:
    cell_types = (
        df.groupby("cell_type", observed=True)[ifn_pattern]
        .median()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    conditions = ["ctrl", "stim"]
    colors = {"ctrl": "#6B7280", "stim": "#2B8CBE"}
    positions = []
    labels = []
    data = []
    color_list = []
    for i, cell_type in enumerate(cell_types):
        base = i * 3
        for j, condition in enumerate(conditions):
            vals = df.loc[
                (df["cell_type"] == cell_type) & (df["condition"] == condition),
                ifn_pattern,
            ].dropna()
            positions.append(base + j)
            labels.append(cell_type if j == 0 else "")
            data.append(vals.to_numpy())
            color_list.append(colors[condition])

    fig, ax = plt.subplots(figsize=(12, 5.8))
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.7,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
    )
    for patch, color in zip(bp["boxes"], color_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor("#111827")

    ax.set_xticks([i * 3 + 0.5 for i in range(len(cell_types))])
    ax.set_xticklabels(cell_types, rotation=35, ha="right")
    ax.set_ylabel(f"{ifn_pattern} activity")
    ax.set_title("Selected R K6 IFN-associated pattern activity by cell type and condition")
    handles = [
        plt.Line2D([0], [0], color=colors["ctrl"], lw=8, label="ctrl"),
        plt.Line2D([0], [0], color=colors["stim"], lw=8, label="stim"),
    ]
    ax.legend(handles=handles, title="condition", frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    save_figure(fig, "ifn_pattern_by_celltype_split_condition.png", dpi=180)


def save_ifn_by_celltype_all(df: pd.DataFrame, ifn_pattern: str) -> None:
    cell_types = (
        df.groupby("cell_type", observed=True)[ifn_pattern]
        .median()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    data = [
        df.loc[df["cell_type"] == cell_type, ifn_pattern].dropna().to_numpy()
        for cell_type in cell_types
    ]
    fig, ax = plt.subplots(figsize=(11, 5.4))
    bp = ax.boxplot(
        data,
        widths=0.65,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#2B8CBE")
        patch.set_alpha(0.72)
        patch.set_edgecolor("#111827")
    ax.set_xticklabels(cell_types, rotation=35, ha="right")
    ax.set_ylabel(f"{ifn_pattern} activity")
    ax.set_title("Selected R K6 IFN-associated pattern activity by cell type")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    save_figure(fig, "ifn_pattern_by_celltype_all.png", dpi=180)


def save_mean_pattern_by_condition(df: pd.DataFrame) -> None:
    patterns = pattern_columns(df)
    condition_order = ["ctrl", "stim"]
    means = (
        df.groupby("condition", observed=True)[patterns]
        .mean()
        .reindex(condition_order)
    )
    values = means.to_numpy()
    fig, ax = plt.subplots(figsize=(9, 3.8))
    im = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(patterns)))
    ax.set_xticklabels(patterns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(condition_order)))
    ax.set_yticklabels(condition_order)
    ax.set_title("Selected R K6 mean pattern activity by condition")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Mean activity")
    fig.tight_layout()
    save_figure(fig, "mean_pattern_by_condition.png", dpi=180)


def main() -> None:
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    if LEGACY_IMG_DIR is not None:
        LEGACY_IMG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CELL_ACTIVITIES)
    ifn_pattern = selected_ifn_pattern()
    save_pattern_condition_overview(ifn_pattern)
    save_umap_pattern_pair_figures(df, ifn_pattern)
    save_ifn_replicate_pairs(ifn_pattern)
    save_ifn_celltype_condition(ifn_pattern)
    save_pattern_celltype_shift_heatmap()
    save_top_gene_directionality(ifn_pattern)
    identity_pattern = selected_identity_pattern(ifn_pattern)
    save_identity_pattern_directionality(identity_pattern)
    save_ifn_by_celltype_split(df, ifn_pattern)
    save_ifn_by_celltype_all(df, ifn_pattern)
    save_mean_pattern_by_condition(df)
    if DIRECTIONALITY_HEATMAP.exists():
        target = IMG_DIR / "pattern_gene_directionality_heatmap.png"
        target.write_bytes(DIRECTIONALITY_HEATMAP.read_bytes())
        if LEGACY_IMG_DIR is not None:
            legacy_target = LEGACY_IMG_DIR / "pattern_gene_directionality_heatmap.png"
            legacy_target.write_bytes(DIRECTIONALITY_HEATMAP.read_bytes())
    print(f"Selected K6 IFN-associated pattern: {ifn_pattern}")
    print(f"Wrote figures to {IMG_DIR}")
    print(f"Wrote source tables to {SOURCE_TABLE_DIR}")


if __name__ == "__main__":
    main()
