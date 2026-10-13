#!/usr/bin/env python3
"""Regression tests for source-wide CPM; no CoGAPS fit or production writes."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import directionality_normalization as norm
import pattern_directionality_analysis as direction


class SourceCountsTests(unittest.TestCase):
    def setUp(self):
        obs = pd.DataFrame({"replicate": ["d1"] * 4, "label": ["ctrl", "ctrl", "stim", "stim"],
                            "cell_type": ["B"] * 4}, index=["a", "b", "c", "d"])
        self.source = ad.AnnData(sparse.csr_matrix([[1, 2, 3], [2, 5, 7], [2, 1, 10], [4, 3, 20]]),
                                obs=obs, var=pd.DataFrame(index=["g-1", "g2", "g3"]))
        self.prepared = self.source[:, :2].copy()
        self.prepared.obs["condition"] = self.prepared.obs.label
        self.prepared.layers["counts"] = self.prepared.X.copy()
        self.prepared.X = self.prepared.X.astype(float) * 0.5

    def totals(self):
        table = norm.verify_source_counts(self.source, self.prepared)
        return table.set_index("cell_id").source_total_counts

    def test_hand_calculation_and_identifier_reordering(self):
        totals = self.totals().iloc[::-1]
        pb = direction.make_pseudobulk(self.prepared, ["g-1"], ["replicate", "condition"], 2, totals)
        pb = direction.add_cpm(pb, ["g-1"])
        np.testing.assert_array_equal(pb.library_size_total_counts, [20, 40])
        np.testing.assert_array_equal(pb["g-1"], [3, 6])
        np.testing.assert_array_equal(pb["g-1__cpm"], [150000, 150000])
        pairs, _ = direction.paired_directionality(pb, ["g-1"], ["replicate"], {})
        self.assertEqual(pairs.log2fc_stim_vs_ctrl.iloc[0], 0)
        reordered = norm.verify_source_counts(self.source[::-1, ::-1].copy(), self.prepared)
        np.testing.assert_array_equal(reordered.source_total_counts, [6, 14, 13, 27])

    def test_no_transformed_fallback(self):
        del self.prepared.layers["counts"]
        with self.assertRaisesRegex(ValueError, "fallback"):
            direction.get_counts_matrix(self.prepared)

    def test_invalid_counts(self):
        for value in [-1, 0.5, np.nan, np.inf, float(2**53)]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                norm.raw_counts([[value]])

    def test_duplicate_and_missing_identifiers(self):
        for ids in [["a", "a"], ["a", ""], ["a", None]]:
            with self.subTest(ids=ids), self.assertRaises(ValueError):
                norm.require_ids(ids, "fixture")
        self.source.obs_names = ["a", "b", "c", "extra"]
        with self.assertRaisesRegex(ValueError, "Cell sets"):
            self.totals()

    def test_changed_numerator_is_rejected(self):
        self.prepared.layers["counts"][0, 0] = 9
        with self.assertRaisesRegex(ValueError, "disagree"):
            self.totals()

    def test_changed_metadata_is_rejected(self):
        self.prepared.obs.loc["a", "replicate"] = "wrong"
        with self.assertRaisesRegex(ValueError, "metadata mismatch"):
            self.totals()

    def test_missing_grouping_metadata_is_rejected(self):
        self.prepared.obs.loc["a", "condition"] = None
        with self.assertRaisesRegex(ValueError, "Missing grouping"):
            self.totals()

    def test_zero_source_library_is_rejected(self):
        self.source.X[0, :] = 0
        self.prepared.layers["counts"][0, :] = 0
        with self.assertRaisesRegex(ValueError, "positive"):
            self.totals()

    def test_invalid_denominators_are_not_zero_cpm(self):
        for value in [0, -1, np.nan, np.inf]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                direction.add_cpm(pd.DataFrame({"library_size_total_counts": [value], "g": [1]}), ["g"])

    def test_pair_duplicate_and_nonfinite_rejected(self):
        pb = pd.DataFrame({"replicate": ["d1", "d1"], "condition": ["ctrl", "ctrl"],
                           "g__cpm": [1, 2], "n_cells": [20, 20]})
        with self.assertRaisesRegex(ValueError, "at most one"):
            direction.paired_directionality(pb, ["g"], ["replicate"], {})
        pb.loc[1, "condition"] = "stim"
        pb.loc[1, "g__cpm"] = np.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            direction.paired_directionality(pb, ["g"], ["replicate"], {})

    def test_totals_must_cover_same_cells(self):
        with self.assertRaisesRegex(ValueError, "exactly"):
            direction.make_pseudobulk(self.prepared, ["g-1"], ["replicate", "condition"], 2, self.totals().iloc[:3])

    def test_cell_threshold_is_explicit(self):
        obs = self.prepared.obs.copy()
        obs.loc["d", "cell_type"] = "T"
        coverage = norm.group_eligibility(obs, ["replicate", "condition", "cell_type"], 2)
        self.assertEqual(coverage.paired_eligible.sum(), 0)
        self.assertEqual(sum(coverage.exclusion_reason == "below_min_cells"), 2)
        self.assertEqual(sum(coverage.exclusion_reason == "counterpart_missing_or_below_min_cells"), 1)
        with self.assertRaisesRegex(ValueError, "No pseudobulk"):
            direction.make_pseudobulk(self.prepared, ["g-1"], ["replicate", "condition"], 3, self.totals())

    def test_manifest_integrity_and_prepared_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared_path = root / "prepared.h5ad"
            self.prepared.write_h5ad(prepared_path)
            norm.verify_source_counts(self.source, self.prepared).to_csv(root / "totals.csv", index=False)
            pd.DataFrame({"gene": self.source.var_names}).to_csv(root / "source.csv", index=False)
            pd.DataFrame({"gene": self.prepared.var_names}).to_csv(root / "modeled.csv", index=False)
            manifest = {"schema_version": norm.SCHEMA, "retained_counts_match_source": True,
                        "prepared": norm.file_record(prepared_path), "n_cells": 4,
                        "source_gene_count": 3, "modeled_gene_count": 2}
            for key, name in (("totals", "totals.csv"), ("source_genes", "source.csv"), ("modeled_genes", "modeled.csv")):
                manifest[key] = norm.file_record(root / name)
            path = root / "manifest.json"
            path.write_text(json.dumps(manifest))
            actual, _ = norm.load_normalization(path, prepared_path, self.prepared)
            np.testing.assert_array_equal(actual, [6, 14, 13, 27])
            with (root / "totals.csv").open("a") as stream:
                stream.write("\n")
            with self.assertRaisesRegex(ValueError, "Checksum"):
                norm.load_normalization(path, prepared_path, self.prepared)
            with prepared_path.open("ab") as stream:
                stream.write(b"changed")
            with self.assertRaisesRegex(ValueError, "Checksum"):
                norm.load_normalization(path, prepared_path, self.prepared)

    def test_output_overwrite_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "out"
            norm.fresh_output(path)
            (path / "keep.txt").write_text("keep")
            with self.assertRaises(FileExistsError):
                norm.fresh_output(path)
            self.assertEqual((path / "keep.txt").read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
