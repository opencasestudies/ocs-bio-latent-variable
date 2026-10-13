#!/usr/bin/env python3
"""Check saved directionality imports without source data, full models, or fitting."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


BUNDLE_RELATIVE = Path("data/processed/selected_model_k6/directionality_source_wide")
METHODS = {"hvg_cpm", "source_cpm", "tmm_cpm"}
LEGACY_FIELDS = {"evidence_strength", "sign_consistency", "wilcoxon_pvalue", "wilcoxon_fdr"}
REPORT_FIELDS = {
    "pattern", "gene", "rank", "cogaps_weight", "scope", "cell_type",
    "normalization", "mean_log2_contrast", "median_log2_contrast", "n_pairs",
    "n_positive", "n_negative", "n_zero", "n_agree_mean", "agreement_with_mean",
    "mean_direction", "donor_support", "eligible_donors", "normalization_status",
    "hvg_cpm", "source_cpm", "tmm_cpm", "direction_stable",
}


def contained(root, relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"Expected a bundle-relative path: {relative}")
    resolved = (root / path).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f"Path leaves the bundle: {relative}")
    return resolved


def verify(path, record):
    if not path.is_file() or path.stat().st_size != record["size_bytes"]:
        raise ValueError(f"Missing file or byte-size mismatch: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != record["sha256"]:
        raise ValueError(f"SHA-256 mismatch: {path}")


def read_table(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames or []
        rows = list(reader)
    if not columns or len(columns) != len(set(columns)) or not rows:
        raise ValueError(f"Empty table or duplicate columns: {path}")
    if LEGACY_FIELDS.intersection(columns):
        raise ValueError(f"Legacy confidence/test columns: {path}")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Malformed CSV row: {path}")
    return columns, rows


def sign(value):
    return 1 if value > 0 else -1 if value < 0 else 0


def validate_report(path):
    columns, rows = read_table(path)
    if not REPORT_FIELDS.issubset(columns):
        raise ValueError(f"Missing donor/sensitivity fields: {path}")
    keys = set()
    for row in rows:
        key = tuple(row[name] for name in ("pattern", "gene", "scope", "cell_type"))
        if key in keys:
            raise ValueError(f"Duplicate report key: {key}")
        keys.add(key)
        n, pos, neg, zero = (int(row[name]) for name in ("n_pairs", "n_positive", "n_negative", "n_zero"))
        donors = row["eligible_donors"].split(";") if row["eligible_donors"] else []
        if min(n, pos, neg, zero) < 0 or pos + neg + zero != n or len(donors) != n or len(set(donors)) != n:
            raise ValueError(f"Inconsistent donor counts: {key}")
        if row["normalization"] != "source_cpm":
            raise ValueError(f"Primary report must use source-wide CPM: {key}")
        if not n:
            if row["donor_support"] != "not estimable" or row["normalization_status"] != "not estimable":
                raise ValueError(f"Missing not-estimable label: {key}")
            continue
        mean, median = (float(row[name]) for name in ("mean_log2_contrast", "median_log2_contrast"))
        method_means = [float(row[name]) for name in sorted(METHODS)]
        if not all(math.isfinite(value) for value in (mean, median, *method_means)):
            raise ValueError(f"Nonfinite estimate: {key}")
        if not math.isclose(mean, float(row["source_cpm"]), abs_tol=1e-12):
            raise ValueError(f"Primary mean differs from source CPM: {key}")
        stable = len({sign(value) for value in method_means}) == 1
        label = "same mean sign across methods" if stable else "normalization-sensitive mean sign"
        if row["direction_stable"].lower() != str(stable).lower() or row["normalization_status"] != label:
            raise ValueError(f"Incorrect normalization flag: {key}")
        direction = "positive" if mean > 0 else "negative" if mean < 0 else "zero"
        if row["mean_direction"] != direction:
            raise ValueError(f"Incorrect mean direction: {key}")
        expected_support = (
            "single donor pair; no across-donor replication" if n == 1 else
            "mixed donor directions" if pos and neg else
            "no observed contrast at this scale" if zero == n else
            f"{pos or neg} of {n} {'positive' if pos else 'negative'}"
        )
        if row["donor_support"] != expected_support:
            raise ValueError(f"Incorrect donor description: {key}")
        if mean:
            agree = pos if mean > 0 else neg
            if float(row["n_agree_mean"]) != agree or not math.isclose(float(row["agreement_with_mean"]), agree / n, abs_tol=1e-12):
                raise ValueError(f"Incorrect agreement with mean: {key}")
        elif row["agreement_with_mean"] not in ("", "NA", "NaN", "nan"):
            raise ValueError(f"Zero-mean agreement must be undefined: {key}")
    return len(rows)


def check_bundle(bundle, prepared=None):
    bundle = Path(bundle).resolve()
    manifest = json.loads((bundle / "bundle_manifest.json").read_text())
    if manifest["schema_version"] != "directionality-bundle-v1" or manifest["reporting_version"] != "donor-support-v1":
        raise ValueError("Unsupported directionality bundle schema")
    indexed = {record["path"]: record for record in manifest["files"]}
    if len(indexed) != len(manifest["files"]):
        raise ValueError("Duplicate bundle manifest entries")
    required_manifests = {
        "normalization/normalization_manifest.json", "r/directionality_method.json",
        "python/directionality_method.json", "sensitivity/sensitivity_method.json",
    }
    if not required_manifests.issubset(indexed):
        raise ValueError("Missing indexed component manifest")
    for name, record in indexed.items():
        path = contained(bundle, name)
        verify(path, record)
        if path.suffix == ".csv":
            columns, rows = read_table(path)
            if columns != record["columns"] or len(rows) != record["rows"]:
                raise ValueError(f"CSV schema/row-count mismatch: {name}")
    norm = json.loads((bundle / "normalization/normalization_manifest.json").read_text())
    if (norm["source_gene_count"], norm["modeled_gene_count"], norm["n_cells"]) != (15706, 3000, 24673):
        raise ValueError("Bundle does not describe the selected prepared dataset")
    for field in ("totals", "source_genes", "modeled_genes"):
        record = norm[field]
        name = f"normalization/{record['filename']}"
        if name not in indexed:
            raise ValueError(f"Normalization input not indexed: {name}")
        verify(contained(bundle, name), record)
    if prepared is not None:
        verify(Path(prepared), norm["prepared"])
    components = {"r": "directionality_method.json", "python": "directionality_method.json", "sensitivity": "sensitivity_method.json"}
    for directory, filename in components.items():
        method = json.loads((bundle / directory / filename).read_text())
        if method["normalization"] != norm or method["reporting_version"] != "donor-support-v1":
            raise ValueError(f"Mismatched normalization/reporting contract: {directory}")
        for record in method["outputs"]:
            name = f"{directory}/{record['filename']}"
            if name not in indexed:
                raise ValueError(f"Output not indexed: {name}")
            verify(contained(bundle, name), record)
    sensitivity = json.loads((bundle / "sensitivity/sensitivity_method.json").read_text())
    if sensitivity["schema_version"] != "directionality-sensitivity-v2" or set(sensitivity["methods"]) != METHODS or sensitivity["pseudocount_CPM"] != 0.5:
        raise ValueError("Unsupported sensitivity calculation")
    reports = {}
    for language in ("r", "python"):
        for scope in ("global", "by_celltype", "identity"):
            name = f"sensitivity/{language}_gene_report_{scope}.csv"
            if name not in indexed:
                raise ValueError(f"Missing canonical report: {name}")
            reports[name] = validate_report(bundle / name)
    return {
        "ready": True, "integration_state": manifest["integration_state"],
        "files_checked": len(indexed), "canonical_report_rows": reports,
        "prepared_identity_checked": prepared is not None,
        "scope": "Saved-file integrity and reporting contract; not a model fit, scientific release approval, or completed lesson integration.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path(__file__).resolve().parents[2] / BUNDLE_RELATIVE)
    parser.add_argument("--prepared", type=Path, help="Optionally verify the included prepared H5AD identity.")
    args = parser.parse_args()
    try:
        result = check_bundle(args.bundle, args.prepared)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
