#!/usr/bin/env python3
"""Recreate the prepared PBMC dataset from a larger source H5AD file.

This script is for the optional full reproduction path. The main learner path
uses the prepared H5AD and lightweight CoGAPS outputs already included with the
case study.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-h5ad",
        type=Path,
        default=Path("data/source/kang_counts_25k.h5ad"),
        help="Larger source counts H5AD file.",
    )
    parser.add_argument(
        "--out-h5ad",
        type=Path,
        default=Path("data/processed/input/preprocessed_cells_hvg3000.h5ad"),
        help="Prepared cells x genes H5AD file.",
    )
    parser.add_argument(
        "--cogaps-input-h5ad",
        type=Path,
        default=Path("data/processed/input/cogaps_input_genesxcells_hvg3000_float64.h5ad"),
        help="Optional dense genes x cells PyCoGAPS input.",
    )
    parser.add_argument(
        "--config-json",
        type=Path,
        default=Path("data/processed/input/preprocess_config_hvg3000.json"),
        help="Preprocessing settings record.",
    )
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--target-sum", type=float, default=10000.0)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--hvg-flavor", default="seurat_v3")
    return parser.parse_args()


def build_cogaps_input(preprocessed_cells: ad.AnnData) -> ad.AnnData:
    """Convert cells x genes AnnData into a dense genes x cells file."""
    cogaps_input = preprocessed_cells.T.copy()
    X = cogaps_input.X.toarray() if sparse.issparse(cogaps_input.X) else cogaps_input.X
    cogaps_input.X = np.asarray(X, dtype=np.float64)
    cogaps_input.var = preprocessed_cells.obs.copy()
    cogaps_input.var_names = preprocessed_cells.obs_names.astype(str).copy()
    cogaps_input.obs = preprocessed_cells.var.copy()
    cogaps_input.obs_names = preprocessed_cells.var_names.astype(str).copy()
    return cogaps_input


def main() -> None:
    args = parse_args()

    try:
        import scanpy as sc
    except ImportError as exc:
        raise SystemExit(
            "This full-reproduction preprocessing script requires scanpy. "
            "Use the validated case-study image or install scanpy in your environment."
        ) from exc

    if not args.source_h5ad.exists():
        raise FileNotFoundError(
            f"Source H5AD not found: {args.source_h5ad}. "
            "Download or mount the larger source file before running this script."
        )

    adata = ad.read_h5ad(args.source_h5ad)

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        flavor=args.hvg_flavor,
        layer="counts",
        subset=False,
    )

    if "highly_variable" not in adata.var:
        raise RuntimeError("HVG selection did not create adata.var['highly_variable'].")

    adata_hvg = adata[:, adata.var["highly_variable"].to_numpy()].copy()

    args.out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    args.cogaps_input_h5ad.parent.mkdir(parents=True, exist_ok=True)
    args.config_json.parent.mkdir(parents=True, exist_ok=True)

    adata_hvg.write_h5ad(args.out_h5ad)

    cogaps_input = build_cogaps_input(adata_hvg)
    cogaps_input.write_h5ad(args.cogaps_input_h5ad)

    config = {
        "raw_h5ad": args.source_h5ad.name,
        "min_cells": args.min_cells,
        "target_sum": args.target_sum,
        "hvg_flavor": args.hvg_flavor,
        "n_top_genes": args.n_top_genes,
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "prepared_h5ad": str(args.out_h5ad),
        "cogaps_input_h5ad": str(args.cogaps_input_h5ad),
    }
    args.config_json.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote prepared PBMC dataset: {args.out_h5ad}")
    print(f"Wrote CoGAPS input: {args.cogaps_input_h5ad}")
    print(f"Wrote preprocessing record: {args.config_json}")


if __name__ == "__main__":
    main()
