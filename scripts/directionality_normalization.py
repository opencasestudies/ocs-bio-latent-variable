#!/usr/bin/env python3
"""Build and verify barcode-keyed source-wide library totals for directionality.

This is a reproduction step, not a model fit. The archived source stores counts
in X; the prepared object's X is transformed and must never substitute for counts.
"""

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


SCHEMA = "source-wide-counts-v1"
ROOT = Path(__file__).resolve().parents[1]


def file_record(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return {"filename": path.name, "size_bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def verify_file(path, record):
    actual = file_record(path)
    if any(actual[key] != record[key] for key in ("size_bytes", "sha256")):
        raise ValueError(f"Checksum/size mismatch: {path}")


def fresh_output(path):
    path = Path(path)
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise FileExistsError(f"Use an empty output directory; existing files are protected: {path}")
    path.mkdir(parents=True, exist_ok=True)


def require_ids(values, label):
    values = pd.Index(values)
    if values.hasnans or values.has_duplicates or any(not str(x).strip() for x in values):
        raise ValueError(f"{label}: identifiers must be unique and nonmissing")
    return values.astype(str)


def raw_counts(matrix):
    counts = sparse.csr_matrix(matrix, dtype=np.float64)
    values = counts.data
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Raw counts must be finite, nonnegative integers")
    # R uses doubles; stay within exact-integer range in both implementations.
    if values.sum() >= 2**53:
        raise ValueError("Count sum exceeds the exact-integer range")
    return counts.astype(np.int64)


def prepared_counts(prepared):
    if "counts" not in prepared.layers:
        raise ValueError("Prepared layers/counts is required; transformed X is not a raw-count fallback")
    return raw_counts(prepared.layers["counts"])


def metadata(prepared):
    columns = ["replicate", "condition", "cell_type"]
    obs = prepared.obs[columns].copy()
    for column in columns:
        if obs[column].isna().any() or obs[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing grouping metadata: {column}")
        obs[column] = obs[column].astype(str)
    if set(obs.condition) != {"ctrl", "stim"}:
        raise ValueError("Expected exactly ctrl and stim conditions")
    return obs


def verify_source_counts(source, prepared):
    source_cells = require_ids(source.obs_names, "source cells")
    cells = require_ids(prepared.obs_names, "prepared cells")
    source_genes = require_ids(source.var_names, "source genes")
    genes = require_ids(prepared.var_names, "prepared genes")
    if set(source_cells) != set(cells) or not set(genes).issubset(source_genes):
        raise ValueError("Cell sets must match exactly and all modeled genes must occur in the source")
    cell_index = source_cells.get_indexer(cells)
    gene_index = source_genes.get_indexer(genes)
    counts = prepared_counts(prepared)
    full_counts = raw_counts(source.X)[cell_index]
    if (full_counts[:, gene_index] != counts).nnz:
        raise ValueError("Prepared counts disagree with the source at matched cell/gene identifiers")
    obs = metadata(prepared)
    for target, origin in (("replicate", "replicate"), ("condition", "label"), ("cell_type", "cell_type")):
        values = source.obs.iloc[cell_index][origin]
        if values.isna().any() or not np.array_equal(values.astype(str), obs[target]):
            raise ValueError(f"Source/prepared metadata mismatch: {target}")
    totals = np.asarray(full_counts.sum(axis=1)).ravel()
    if (totals <= 0).any():
        raise ValueError("Source library totals must be positive for every retained cell")
    table = obs.reset_index(drop=True)
    table.insert(0, "cell_id", cells)
    table["source_total_counts"] = totals
    table["modeled_total_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    return table


def load_normalization(manifest_path, prepared_path, prepared):
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if manifest["schema_version"] != SCHEMA or manifest["retained_counts_match_source"] is not True:
        raise ValueError("Unsupported or unverified normalization manifest")
    verify_file(prepared_path, manifest["prepared"])
    tables = {}
    for key in ("totals", "source_genes", "modeled_genes"):
        path = manifest_path.parent / manifest[key]["filename"]
        verify_file(path, manifest[key])
        tables[key] = pd.read_csv(path, dtype=str, keep_default_na=False)
    cells = require_ids(prepared.obs_names, "prepared cells")
    genes = require_ids(prepared.var_names, "prepared genes")
    table = tables["totals"]
    ids = require_ids(table.cell_id, "normalization cells")
    if set(ids) != set(cells) or len(cells) != manifest["n_cells"]:
        raise ValueError("Normalization cell set differs from prepared data")
    table = table.set_index("cell_id").loc[cells]
    source_genes = require_ids(tables["source_genes"].gene, "normalization source genes")
    modeled_genes = require_ids(tables["modeled_genes"].gene, "normalization modeled genes")
    if (set(modeled_genes) != set(genes) or not set(genes).issubset(source_genes)
            or len(source_genes) != manifest["source_gene_count"]
            or len(genes) != manifest["modeled_gene_count"]):
        raise ValueError("Normalization gene universes disagree")
    counts = prepared_counts(prepared)
    expected_hvg = np.asarray(counts.sum(axis=1)).ravel()
    if not np.array_equal(pd.to_numeric(table.modeled_total_counts), expected_hvg):
        raise ValueError("Normalization modeled totals disagree with prepared counts")
    obs = metadata(prepared)
    if not np.array_equal(table[obs.columns].to_numpy(), obs.to_numpy()):
        raise ValueError("Normalization grouping metadata disagree")
    totals = pd.to_numeric(table.source_total_counts).to_numpy(dtype=np.float64)
    if (not np.isfinite(totals).all() or (totals <= 0).any()
            or (totals != np.floor(totals)).any() or (totals < expected_hvg).any()
            or totals.sum() >= 2**53):
        raise ValueError("Invalid source-wide library totals")
    return pd.Series(totals, index=cells), manifest


def group_eligibility(obs, group_cols, min_cells):
    groups = obs.groupby(group_cols, observed=True, dropna=False).size().rename("n_cells").reset_index()
    groups["passes_min_cells"] = groups.n_cells >= min_cells
    pair_cols = [c for c in group_cols if c != "condition"]
    eligible = groups.groupby(pair_cols, observed=True)["passes_min_cells"].transform("sum")
    groups["paired_eligible"] = groups.passes_min_cells & eligible.eq(2)
    groups["exclusion_reason"] = np.where(
        ~groups.passes_min_cells, "below_min_cells",
        np.where(groups.paired_eligible, "included", "counterpart_missing_or_below_min_cells"),
    )
    return groups


def write_method_manifest(out, normalization, model_path, script_path, min_cells, pseudocount):
    out = Path(out)
    document = {
        "normalization": normalization,
        "model": file_record(model_path), "script": file_record(script_path),
        "normalization_helper": file_record(__file__),
        "numerator": "Sum of prepared layers/counts for one gene in the exact grouped cells",
        "denominator": "Sum of source-wide per-cell raw totals for those same cells",
        "grouping_global": ["replicate", "condition"],
        "grouping_cell_type": ["replicate", "condition", "cell_type"],
        "min_cells_per_pseudobulk": min_cells,
        "pseudocount": pseudocount, "pseudocount_unit": "CPM",
        "contrast": "log2(stim_CPM + 0.5) - log2(ctrl_CPM + 0.5)",
        "donor_weighting": "Equal weight among eligible donor-paired contrasts",
        "status": "Candidate outputs; not approved for replacing lesson results",
        "reporting_version": "donor-support-v1",
        "reporting_helper": file_record(Path(__file__).with_name("directionality_reporting.py")),
        "caution": "Descriptive donor counts and mean-sign agreement only; no evidence-strength category or Wilcoxon/FDR. Source-only exports mark sensitivity not assessed until the sensitivity step.",
        "outputs": [file_record(p) for p in sorted(out.iterdir()) if p.is_file()],
    }
    (out / "directionality_method.json").write_text(json.dumps(document, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-h5ad", type=Path, required=True)
    parser.add_argument("--preprocessed-h5ad", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--archive-manifest", type=Path, default=ROOT / "data/external/reproduction_files.json")
    args = parser.parse_args()
    archive = json.loads(args.archive_manifest.read_text())
    source_record = next(row for row in archive["files"] if row["id"] == "source")
    verify_file(args.source_h5ad, source_record)
    source = ad.read_h5ad(args.source_h5ad)
    prepared = ad.read_h5ad(args.preprocessed_h5ad)
    totals = verify_source_counts(source, prepared)
    fresh_output(args.outdir)
    totals.to_csv(args.outdir / "source_library_totals.csv", index=False)
    pd.DataFrame({"gene": source.var_names}).to_csv(args.outdir / "source_genes.csv", index=False)
    pd.DataFrame({"gene": prepared.var_names}).to_csv(args.outdir / "modeled_genes.csv", index=False)
    manifest = {
        "schema_version": SCHEMA, "source": file_record(args.source_h5ad),
        "prepared": file_record(args.preprocessed_h5ad), "exporter": file_record(__file__),
        "archive_manifest": file_record(args.archive_manifest),
        "archive_record_url": archive["record_url"],
        "source_count_location": "X", "prepared_count_location": "layers/counts",
        "source_gene_count": source.n_vars, "modeled_gene_count": prepared.n_vars,
        "n_cells": prepared.n_obs, "retained_counts_match_source": True,
        "normalization_universe": "All genes in the supplied source object; not the entire transcriptome",
        "totals": file_record(args.outdir / "source_library_totals.csv"),
        "source_genes": file_record(args.outdir / "source_genes.csv"),
        "modeled_genes": file_record(args.outdir / "modeled_genes.csv"),
    }
    (args.outdir / "normalization_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Verified {prepared.n_obs} cells; {prepared.n_vars} modeled genes; {source.n_vars} source genes.")


if __name__ == "__main__":
    main()
