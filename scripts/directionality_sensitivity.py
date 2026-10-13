#!/usr/bin/env python3
"""Compare matched descriptive contrasts under HVG, source-wide, and TMM CPM.

Identity contrasts pool all non-focal raw counts within donor and condition.
No CoGAPS fit, differential-expression test, or lesson-output replacement occurs.
"""

import argparse
import gzip
import json
from pathlib import Path
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy import sparse

from directionality_normalization import (
    file_record, fresh_output, load_normalization, metadata, raw_counts,
    verify_file, verify_source_counts,
)
from directionality_reporting import REPORTING_VERSION, summarize_donors


METHODS = ("hvg_cpm", "source_cpm", "tmm_cpm")
PSEUDOCOUNT = 0.5
MIN_CELLS = 20
ROOT = Path(__file__).resolve().parents[1]


def build_pseudobulks(counts, obs, modeled_indices, min_cells=MIN_CELLS):
    """Use identical complete donor sets across methods; retain exclusion records."""
    donors = sorted(obs.replicate.unique())
    cell_types = sorted(obs.cell_type.unique())
    cache = {}
    for donor in donors:
        for condition in ("ctrl", "stim"):
            mask = (obs.replicate == donor) & (obs.condition == condition)
            rows = np.flatnonzero(mask.to_numpy())
            cache[donor, condition, "all"] = (np.asarray(counts[rows].sum(axis=0)).ravel(), len(rows))
            for cell_type in cell_types:
                rows = np.flatnonzero((mask & (obs.cell_type == cell_type)).to_numpy())
                cache[donor, condition, cell_type] = (np.asarray(counts[rows].sum(axis=0)).ravel(), len(rows))
    vectors, samples, coverage = [], [], []
    cohorts = [("stimulation_global", "all")]
    cohorts += [(scope, cell_type) for scope in ("stimulation_celltype", "identity") for cell_type in cell_types]
    for cohort_i, (scope, cell_type) in enumerate(cohorts):
        for donor in donors:
            entries = []
            for condition in ("ctrl", "stim"):
                focal, n_focal = cache[donor, condition, cell_type]
                entries.append((condition, "focal", focal, n_focal))
                if scope == "identity":
                    all_counts, n_all = cache[donor, condition, "all"]
                    entries.append((condition, "other", all_counts - focal, n_all - n_focal))
            eligible = all(n >= min_cells and vector.sum() > 0 for _, _, vector, n in entries)
            for condition, role, vector, n_cells in entries:
                if (vector < 0).any():
                    raise ValueError("Pooled-other counts cannot be negative")
                row = {"cohort": f"cohort_{cohort_i:02d}", "scope": scope, "cell_type": cell_type,
                       "replicate": donor, "condition": condition, "role": role, "n_cells": n_cells,
                       "source_library_size": int(vector.sum()),
                       "hvg_library_size": int(vector[modeled_indices].sum()), "pair_eligible": eligible,
                       "exclusion_reason": "included" if eligible else (
                           "below_min_cells" if n_cells < min_cells else "required_partner_below_min_cells")}
                coverage.append(row)
                if eligible:
                    row = {**row, "sample_id": f"sample_{len(samples):04d}"}
                    samples.append(row)
                    vectors.append(vector)
    if not vectors:
        raise ValueError("No complete donor strata pass eligibility")
    return np.column_stack(vectors), pd.DataFrame(samples), pd.DataFrame(coverage)


def paired_contrasts(cpm, samples, genes):
    if cpm.shape != (len(genes), len(samples)) or not np.isfinite(cpm).all() or (cpm < 0).any():
        raise ValueError("CPM matrix must be finite, nonnegative, and aligned to samples and genes")
    records = []
    samples = samples.reset_index(drop=True)
    for (scope, cell_type, donor), group in samples.groupby(["scope", "cell_type", "replicate"], sort=True):
        expected = {("ctrl", "focal"), ("stim", "focal")}
        if scope == "identity":
            expected |= {("ctrl", "other"), ("stim", "other")}
        if group.duplicated(["condition", "role"]).any() or set(zip(group.condition, group.role)) != expected:
            raise ValueError("Every contrast requires exactly its complete donor/condition/role set")
        idx = {(r.condition, r.role): i for i, r in group.iterrows()}
        if scope == "identity":
            by_condition = {}
            for condition in ("ctrl", "stim"):
                focal, other = idx[condition, "focal"], idx[condition, "other"]
                by_condition[condition] = np.log2(cpm[:, focal] + PSEUDOCOUNT) - np.log2(cpm[:, other] + PSEUDOCOUNT)
            values = (by_condition["ctrl"] + by_condition["stim"]) / 2
        else:
            values = np.log2(cpm[:, idx["stim", "focal"]] + PSEUDOCOUNT) - np.log2(cpm[:, idx["ctrl", "focal"]] + PSEUDOCOUNT)
            by_condition = {}
        for g, gene in enumerate(genes):
            records.append({"scope": scope, "cell_type": cell_type, "replicate": donor, "gene": gene,
                            "log2_contrast": values[g],
                            **{f"identity_{condition}_log2_contrast": vals[g] for condition, vals in by_condition.items()}})
    return pd.DataFrame(records)


def summarize_pairs(pairs, coverage=None):
    keys = ["scope", "cell_type", "gene", "normalization"]
    if pairs[keys + ["replicate", "log2_contrast"]].isna().any().any():
        raise ValueError("Missing donor contrasts or grouping metadata")
    groups = {key: group for key, group in pairs.groupby(keys, sort=True)}
    expected = list(groups)
    if coverage is not None:
        expected = [(scope, cell_type, gene, method)
                    for scope, cell_type in coverage[["scope", "cell_type"]].drop_duplicates().itertuples(index=False, name=None)
                    for gene in sorted(pairs.gene.unique()) for method in METHODS]
    rows = []
    for key in expected:
        group = groups.get(key, pairs.iloc[:0])
        x = group.log2_contrast.to_numpy()
        if not np.isfinite(x).all() or group.replicate.duplicated().any():
            raise ValueError("Summaries require finite, unique donor contrasts")
        rows.append({**dict(zip(keys, key)), **summarize_donors(x), "n_donors": len(x),
                     "eligible_donors": ";".join(sorted(group.replicate.astype(str)))})
    return pd.DataFrame(rows)


def compare_methods(summary):
    keys = ["scope", "cell_type", "gene"]
    wide = summary.pivot(index=keys, columns="normalization", values="mean_log2_contrast").reset_index()
    if not set(METHODS).issubset(wide.columns) or summary.groupby(keys).normalization.nunique().ne(len(METHODS)).any():
        raise ValueError("Comparison requires every normalization method")
    counts = summary.groupby(keys).n_donors.agg(["min", "max"])
    if not (counts["min"] == counts["max"]).all():
        raise ValueError("Methods must use the same eligible donor set")
    if "eligible_donors" in summary and summary.groupby(keys).eligible_donors.nunique().ne(1).any():
        raise ValueError("Methods must use the same eligible donor set")
    wide = wide.merge(counts["min"].rename("n_donors"), on=keys, validate="one_to_one")
    estimable = wide.n_donors.gt(0)
    if not np.isfinite(wide.loc[estimable, list(METHODS)]).all().all():
        raise ValueError("Estimable comparisons require finite means for every method")
    signs = np.sign(wide[list(METHODS)].to_numpy())
    wide["direction_stable"] = pd.array((signs == signs[:, :1]).all(axis=1), dtype="boolean")
    wide["normalization_status"] = np.where(wide.direction_stable, "same mean sign across methods", "normalization-sensitive mean sign")
    wide.loc[~estimable, "direction_stable"] = pd.NA
    wide.loc[~estimable, "normalization_status"] = "not estimable"
    wide["mean_range"] = wide[list(METHODS)].max(axis=1) - wide[list(METHODS)].min(axis=1)
    return wide


def write_reports(pairs, coverage, tops, outdir):
    """Write descriptive reports without changing donor contrasts or inferring confidence."""
    summary = summarize_pairs(pairs, coverage)
    comparison = compare_methods(summary)
    keys = ["scope", "cell_type", "gene"]
    flags = comparison.drop(columns="n_donors")
    summary = summary.merge(flags, on=keys, how="left", validate="many_to_one")
    summary.to_csv(outdir / "contrast_summary.csv", index=False)
    comparison.to_csv(outdir / "normalization_comparison.csv", index=False)
    canonical = summary[summary.normalization.eq("source_cpm")]
    canonical[canonical.scope.eq("identity")].to_csv(outdir / "identity_source_summary.csv", index=False)
    canonical[canonical.scope.eq("stimulation_global")].to_csv(outdir / "stimulation_source_summary.csv", index=False)
    # Repeat genes by pattern only for display; they are not independent tests.
    for language, top in tops.items():
        top = top[["pattern", "gene", "rank", "cogaps_weight"]]
        for scope, suffix in (("stimulation_global", "global"), ("stimulation_celltype", "by_celltype"), ("identity", "identity")):
            report = top.merge(canonical[canonical.scope.eq(scope)], on="gene", how="left", validate="many_to_many")
            if report.normalization_status.isna().any():
                raise ValueError("A pattern gene lacks an explicit report")
            report.to_csv(outdir / f"{language}_gene_report_{suffix}.csv", index=False)
    donor_counts = canonical[["scope", "cell_type", "n_pairs", "eligible_donors"]].drop_duplicates()
    if donor_counts.duplicated(["scope", "cell_type"]).any():
        raise ValueError("Donor eligibility unexpectedly differs among genes")
    donor_counts["replication_note"] = donor_counts.n_pairs.map(lambda n: "not estimable" if n == 0 else "single donor pair; no across-donor replication" if n == 1 else "descriptive comparison; sample size is not a confidence score")
    donor_counts.to_csv(outdir / "eligible_donors_by_scope.csv", index=False)
    return summary, comparison


def load_candidate_genes(directory, filename, normalization):
    document = json.loads((directory / "directionality_method.json").read_text())
    for key in ("prepared", "source"):
        if document["normalization"][key]["sha256"] != normalization[key]["sha256"]:
            raise ValueError("Directionality input provenance differs from sensitivity inputs")
    for record in document["outputs"]:
        verify_file(directory / record["filename"], record)
    return pd.read_csv(directory / filename)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-h5ad", type=Path, required=True)
    parser.add_argument("--preprocessed-h5ad", type=Path, required=True)
    parser.add_argument("--normalization-manifest", type=Path, required=True)
    parser.add_argument("--r-directionality-dir", type=Path, required=True)
    parser.add_argument("--python-directionality-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    fresh_output(args.outdir)
    prepared = ad.read_h5ad(args.preprocessed_h5ad)
    _, normalization = load_normalization(args.normalization_manifest, args.preprocessed_h5ad, prepared)
    verify_file(args.source_h5ad, normalization["source"])
    source = ad.read_h5ad(args.source_h5ad)
    verify_source_counts(source, prepared)
    counts = raw_counts(source.X)[source.obs_names.get_indexer(prepared.obs_names)]
    obs = metadata(prepared)
    modeled_indices = source.var_names.get_indexer(prepared.var_names)
    tops = {"r": load_candidate_genes(args.r_directionality_dir, "pattern_top_genes_directionality.csv", normalization),
            "python": load_candidate_genes(args.python_directionality_dir, "pattern_top_genes.csv", normalization)}
    genes = sorted(set().union(*(set(df.gene) for df in tops.values())))
    gene_indices = source.var_names.get_indexer(genes)
    if (gene_indices < 0).any():
        raise ValueError("A selected pattern gene is absent from the source")
    bulk, samples, coverage = build_pseudobulks(counts, obs, modeled_indices)
    with gzip.open(args.outdir / "broad_counts.mtx.gz", "wb") as stream:
        mmwrite(stream, sparse.coo_matrix(bulk))
    samples.to_csv(args.outdir / "pseudobulk_samples.csv", index=False)
    coverage.to_csv(args.outdir / "donor_eligibility.csv", index=False)
    pd.DataFrame({"gene": source.var_names}).to_csv(args.outdir / "source_genes.csv", index=False)
    for language, table in tops.items():
        table.to_csv(args.outdir / f"pattern_genes_{language}.csv", index=False)
    r_script = ROOT / "scripts/calculate_tmm_factors.R"
    run = subprocess.run(["Rscript", str(r_script), "--input-dir", str(args.outdir),
                          "--outdir", str(args.outdir / "tmm")], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (args.outdir / "tmm_run.log").write_text(run.stdout)
    if run.returncode:
        raise RuntimeError(f"edgeR TMM failed; see {args.outdir / 'tmm_run.log'}")
    factors = pd.read_csv(args.outdir / "tmm/tmm_factors.csv").set_index("sample_id")
    if not factors.index.is_unique or set(factors.index) != set(samples.sample_id):
        raise ValueError("TMM factors must cover exactly the included pseudobulks")
    factors = factors.loc[samples.sample_id]
    np.testing.assert_array_equal(factors.source_library_size, samples.source_library_size)
    methods = {"hvg_cpm": samples.hvg_library_size.to_numpy(float),
               "source_cpm": samples.source_library_size.to_numpy(float),
               "tmm_cpm": factors.effective_library_size.to_numpy(float)}
    all_pairs = []
    for method, denominator in methods.items():
        if not np.isfinite(denominator).all() or (denominator <= 0).any():
            raise ValueError(f"Invalid denominator for {method}")
        cpm = bulk[gene_indices] / denominator[None, :] * 1e6
        pairs = paired_contrasts(cpm, samples, genes)
        pairs["normalization"] = method
        all_pairs.append(pairs)
    pairs = pd.concat(all_pairs, ignore_index=True)
    pairs.to_csv(args.outdir / "donor_contrasts.csv", index=False)
    summary, comparison = write_reports(pairs, coverage, tops, args.outdir)
    method = {
        "schema_version": "directionality-sensitivity-v2", "reporting_version": REPORTING_VERSION,
        "reporting_helper": file_record(ROOT / "scripts/directionality_reporting.py"), "status": "Candidate; not a lesson release",
        "normalization": normalization, "script": file_record(__file__), "tmm_script": file_record(r_script),
        "methods": list(METHODS), "pseudocount_CPM": PSEUDOCOUNT, "min_cells_per_group": MIN_CELLS,
        "identity": "Focal versus ALL other cells pooled as raw counts within donor/condition; both groups >=20 in both conditions; mean of two condition log contrasts per donor then equal mean across donors",
        "stimulation": "Same complete paired-donor sets for all methods; source or HVG library totals, or source totals multiplied by edgeR TMM factors",
        "tmm_cohorts": "One global cohort; one stimulation cohort per cell type; one focal-versus-other cohort per focal type. Complete eligible donors only. All 15706 source genes supplied; all-zero genes omitted internally by edgeR.",
        "hvg_identity_boundary": "Uses the legacy 3000-gene denominator with REPAIRED pooling/pairing and 0.5 pseudocount; does not reproduce the old figure's estimator",
        "scope_boundary": "No absolute RNA claim, significance test, model refit, or confidence category. Identity TMM compares different cell populations and its reference assumption may be less plausible.",
        "outputs": [file_record(p) for p in sorted(args.outdir.iterdir()) if p.is_file()],
        "tmm": json.loads((args.outdir / "tmm/tmm_method.json").read_text()),
    }
    (args.outdir / "sensitivity_method.json").write_text(json.dumps(method, indent=2) + "\n")
    print(comparison.groupby("scope").direction_stable.agg(["size", "sum"]).to_string())


if __name__ == "__main__":
    main()
