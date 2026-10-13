#!/usr/bin/env python3
"""Verify the publication inventory without downloads, fitting, or file changes."""

import argparse
import json
from pathlib import Path

from check_directionality_bundle import BUNDLE_RELATIVE, contained, verify


FIGURE_MANIFEST = Path("data/processed/figures/directionality_figure_manifest.json")
ARCHIVE_INDEX = Path("data/external/reproduction_files.json")


def check_manifest(root):
    root = Path(root).resolve()
    manifest = json.loads((root / "data/artifact_manifest.json").read_text())
    if manifest["schema_version"] != "1.0" or manifest.get("staged_artifact_bundles"):
        raise ValueError("Expected the current inventory without staged artifact bundles")
    included = {}
    for entry in manifest["included_files"]:
        path = contained(root, entry["path"])
        if path in included:
            raise ValueError(f"Duplicate inventory path: {entry['path']}")
        if entry["size_bytes"] <= 0 or entry["github_required"] is not True:
            raise ValueError(f"Invalid included-file record: {entry['path']}")
        verify(path, entry)
        included[path] = entry

    bundle = root / BUNDLE_RELATIVE
    bundle_manifest = json.loads((bundle / "bundle_manifest.json").read_text())
    bundle_paths = [bundle / "bundle_manifest.json", bundle / "README.md"]
    for entry in bundle_manifest["files"]:
        path = contained(bundle, entry["path"])
        verify(path, entry)
        bundle_paths.append(path)
    if not set(bundle_paths).issubset(included):
        raise ValueError("Inventory omits a directionality bundle file")
    if manifest["case_study_model"]["directionality_result_dir"] != BUNDLE_RELATIVE.as_posix():
        raise ValueError("Incorrect directionality result directory")

    index_record = manifest["external_archive_index"]
    if index_record["path"] != ARCHIVE_INDEX.as_posix():
        raise ValueError("Unexpected archive index")
    verify(root / ARCHIVE_INDEX, index_record)
    archive = json.loads((root / ARCHIVE_INDEX).read_text())
    archive_files = {}
    for entry in archive["files"]:
        path = contained(root, entry["target_path"])
        if path in archive_files:
            raise ValueError("Duplicate archive target")
        archive_files[path] = entry
    optional = set()
    for entry in manifest["optional_external_or_local_files"]:
        path = contained(root, entry["recommended_path"])
        if path in optional or path in included or path not in archive_files:
            raise ValueError(f"Duplicate or unindexed optional file: {path}")
        optional.add(path)
        if entry["github_required"] is not False or "present_locally" in entry:
            raise ValueError("Optional records must not depend on the author's local files")
        expected = {key: archive_files[path][key] for key in ("size_bytes", "sha256", "upload_filename")}
        expected.update({key: archive[key] for key in ("record_url", "doi")})
        if any(entry.get(key) != value for key, value in expected.items()):
            raise ValueError(f"Optional file disagrees with the archive index: {path}")

    figure_path = root / FIGURE_MANIFEST
    if figure_path not in included:
        raise ValueError("Inventory omits the figure manifest")
    figure = json.loads(figure_path.read_text())
    # The original generator is historical provenance, not a current-file claim.
    for entry in [figure["validated_generator"], *figure["validated_against"], *figure["outputs"]]:
        verify(contained(root, entry["path"]), entry)

    legacy_names = {"pattern_direction_summary.csv", "pattern_top_genes_directionality.csv",
                    "pseudobulk_counts_global.csv", "pseudobulk_counts_by_celltype.csv"}
    for folder in ("data/processed/selected_model_k6/r", "data/processed/selected_model_k6/python",
                   "data/results", "data/results_python", "data/results_r"):
        for path in (root / folder).glob("*"):
            if path.name in legacy_names or path.name.startswith("pattern_gene_directionality"):
                raise ValueError(f"Superseded directionality file remains: {path}")
    for folder in ("data/processed/figures", "img"):
        if (root / folder / "pattern_gene_directionality_heatmap.png").exists():
            raise ValueError("Superseded directionality heatmap remains")

    return {"ready": True, "included_files_checked": len(included),
            "bundle_files_indexed": len(bundle_paths), "optional_archive_records_checked": len(optional),
            "scope": "Listed file identities and cross-manifest consistency only; no download, model fit, scientific approval, or complete repository inventory."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    try:
        print(json.dumps(check_manifest(args.root), indent=2))
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.exit(1, f"Artifact inventory check failed: {error}\n")
