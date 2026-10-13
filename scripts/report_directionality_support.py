#!/usr/bin/env python3
"""Re-report checksum-verified saved contrasts without recomputing counts or TMM."""

import argparse
import json
from pathlib import Path
import shutil

import pandas as pd

from directionality_normalization import file_record, fresh_output, verify_file
from directionality_reporting import REPORTING_VERSION
from directionality_sensitivity import write_reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sensitivity-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    source = args.sensitivity_dir
    manifest_path = source / "sensitivity_method.json"
    old = json.loads(manifest_path.read_text())
    if old["schema_version"] not in ("directionality-sensitivity-v1", "directionality-sensitivity-v2") or old["pseudocount_CPM"] != 0.5:
        raise ValueError("Unsupported saved contrast contract")
    for record in old["outputs"]:
        verify_file(source / record["filename"], record)
    fresh_output(args.outdir)
    names = ("donor_contrasts.csv", "donor_eligibility.csv", "pattern_genes_r.csv", "pattern_genes_python.csv")
    for name in names:
        if name not in {r["filename"] for r in old["outputs"]}:
            raise ValueError(f"Unverified reporting input: {name}")
        shutil.copyfile(source / name, args.outdir / name)
    pairs = pd.read_csv(source / "donor_contrasts.csv")
    coverage = pd.read_csv(source / "donor_eligibility.csv")
    tops = {language: pd.read_csv(source / f"pattern_genes_{language}.csv") for language in ("r", "python")}
    _, comparison = write_reports(pairs, coverage, tops, args.outdir)
    method = {key: value for key, value in old.items() if key not in ("script", "reporting_helper", "outputs", "reporting_version", "schema_version")}
    method.update({"schema_version": "directionality-sensitivity-v2", "reporting_version": REPORTING_VERSION,
        "script": file_record(__file__), "parent_contrast_manifest": file_record(manifest_path),
        "reporting_helper": file_record(Path(__file__).with_name("directionality_reporting.py")),
        "sensitivity_script": file_record(Path(__file__).with_name("directionality_sensitivity.py")),
        "reporting": "Exact signs of unrounded saved donor values; zero-mean agreement undefined; no confidence categories or p-values. Normalization flag compares three mean signs, not effect magnitude or significance.",
        "outputs": [file_record(p) for p in sorted(args.outdir.iterdir()) if p.is_file()]})
    (args.outdir / "sensitivity_method.json").write_text(json.dumps(method, indent=2) + "\n")
    print(comparison.groupby("scope").direction_stable.agg(["size", "sum"]).to_string())


if __name__ == "__main__":
    main()
