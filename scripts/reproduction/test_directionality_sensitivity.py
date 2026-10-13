#!/usr/bin/env python3

import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from directionality_sensitivity import build_pseudobulks, paired_contrasts, summarize_pairs, compare_methods


class SensitivityTests(unittest.TestCase):
    def fixture(self):
        rows, values = [], []
        for donor in ("d1", "d2"):
            for condition in ("ctrl", "stim"):
                for cell_type, count in (("A", [2, 3]), ("A", [2, 3]), ("B", [1, 19]), ("C", [10, 0])):
                    rows.append(dict(replicate=donor, condition=condition, cell_type=cell_type))
                    values.append(count)
        return sparse.csr_matrix(values), pd.DataFrame(rows)

    def test_other_pool_includes_small_cell_types(self):
        counts, obs = self.fixture()
        bulk, samples, coverage = build_pseudobulks(counts, obs, [0], min_cells=2)
        other = samples[(samples.scope == "identity") & (samples.cell_type == "A") & (samples.role == "other")]
        np.testing.assert_array_equal(bulk[:, other.index], np.tile([[11], [19]], (1, 4)))
        self.assertTrue(other.n_cells.eq(2).all())
        # The raw pooled CPM is 11/30, not the cell-weighted average of 1/20 and 10/10.
        self.assertNotAlmostEqual(11 / 30, (1 / 20 + 10 / 10) / 2)
        self.assertFalse(coverage[(coverage.scope == "identity") & coverage.cell_type.isin(["B", "C"])].pair_eligible.any())

    def test_incomplete_donor_excluded_in_both_conditions(self):
        counts, obs = self.fixture()
        obs.loc[9, "cell_type"] = "B"
        _, samples, coverage = build_pseudobulks(counts, obs, [0], min_cells=2)
        identity = samples[(samples.scope == "identity") & (samples.cell_type == "A")]
        self.assertEqual(set(identity.replicate), {"d1"})
        bad = coverage[(coverage.scope == "identity") & (coverage.cell_type == "A") & (coverage.replicate == "d2")]
        self.assertEqual(len(bad), 4)
        self.assertFalse(bad.pair_eligible.any())

    def test_raw_partition_and_hvg_total(self):
        counts, obs = self.fixture()
        bulk, samples, _ = build_pseudobulks(counts, obs, [0], min_cells=2)
        np.testing.assert_array_equal(bulk.sum(axis=0), samples.source_library_size)
        np.testing.assert_array_equal(bulk[0], samples.hvg_library_size)
        for (donor, condition), group in samples[samples.scope == "identity"].groupby(["replicate", "condition"]):
            actual = bulk[:, group.index].sum(axis=1)
            expected = np.asarray(counts[np.flatnonzero(((obs.replicate == donor) & (obs.condition == condition)).to_numpy())].sum(axis=0)).ravel()
            np.testing.assert_array_equal(actual, expected)

    def test_equal_condition_then_donor_weighting(self):
        samples, cpm = [], []
        for donor, contrasts in (("d1", [2, 4]), ("d2", [-1, 1])):
            for condition, contrast in zip(("ctrl", "stim"), contrasts):
                for role, value in (("focal", 2**contrast - 0.5), ("other", 0.5)):
                    samples.append(dict(scope="identity", cell_type="A", replicate=donor, condition=condition, role=role))
                    cpm.append(value)
        samples = pd.DataFrame(samples)
        pairs = paired_contrasts(np.array([cpm]), samples, ["gene"])
        np.testing.assert_allclose(pairs.log2_contrast, [3, 0], atol=1e-12)
        pairs["normalization"] = "source_cpm"
        summary = summarize_pairs(pairs)
        self.assertEqual(summary.mean_log2_contrast.iloc[0], 1.5)
        with self.assertRaisesRegex(ValueError, "complete"):
            paired_contrasts(np.array([cpm[:-1]]), samples.iloc[:-1], ["gene"])

    def test_nonfinite_and_duplicate_samples(self):
        samples = pd.DataFrame([dict(scope="stimulation_global", cell_type="all", replicate="d1", condition=c, role="focal") for c in ("ctrl", "stim")])
        with self.assertRaisesRegex(ValueError, "finite"):
            paired_contrasts(np.array([[1, np.nan]]), samples, ["gene"])
        samples.loc[1, "condition"] = "ctrl"
        with self.assertRaisesRegex(ValueError, "complete"):
            paired_contrasts(np.array([[1, 2]]), samples, ["gene"])

    def test_agreement_is_with_mean_not_majority(self):
        pairs = pd.DataFrame({"scope": ["identity"] * 3, "cell_type": ["A"] * 3, "gene": ["g"] * 3,
                              "normalization": ["source_cpm"] * 3, "replicate": ["a", "b", "c"], "log2_contrast": [10, -1, -1]})
        result = summarize_pairs(pairs).iloc[0]
        self.assertEqual(result.n_agree_mean, 1)
        self.assertEqual(result.n_negative, 2)

    def test_sign_stability_is_not_effect_or_confidence(self):
        rows = [dict(scope="identity", cell_type="A", gene="g", normalization=m, n_donors=1, mean_log2_contrast=x)
                for m, x in zip(("hvg_cpm", "source_cpm", "tmm_cpm"), (-0.1, 0.05, 0.03))]
        result = compare_methods(pd.DataFrame(rows)).iloc[0]
        self.assertFalse(result.direction_stable)
        self.assertEqual(result.n_donors, 1)
        rows[-1]["n_donors"] = 2
        with self.assertRaisesRegex(ValueError, "same eligible"):
            compare_methods(pd.DataFrame(rows))


if __name__ == "__main__":
    unittest.main()
