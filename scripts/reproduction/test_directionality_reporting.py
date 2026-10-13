#!/usr/bin/env python3

import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from directionality_reporting import summarize_donors
from directionality_sensitivity import summarize_pairs, compare_methods
from pattern_directionality_analysis import paired_directionality


class ReportingTests(unittest.TestCase):
    def test_no_pairs(self):
        r = summarize_donors([])
        self.assertEqual(r["donor_support"], "not estimable")
        self.assertTrue(np.isnan(r["mean_log2_contrast"]))
        self.assertTrue(np.isnan(r["agreement_with_mean"]))

    def test_one_pair_never_implies_replication(self):
        for x in (0, 0.01, -9, 9):
            with self.subTest(x=x):
                r = summarize_donors([x])
                self.assertEqual(r["donor_support"], "single donor pair; no across-donor replication")
                self.assertEqual(r["n_pairs"], 1)

    def test_all_zero(self):
        r = summarize_donors([0, 0, 0])
        self.assertEqual(r["donor_support"], "no observed contrast at this scale")
        self.assertTrue(np.isnan(r["agreement_with_mean"]))
        self.assertTrue(np.isnan(r["n_agree_mean"]))

    def test_positive_mean_negative_majority(self):
        r = summarize_donors([10, -1, -1])
        self.assertGreater(r["mean_log2_contrast"], 0)
        self.assertLess(r["median_log2_contrast"], 0)
        self.assertEqual(r["agreement_with_mean"], 1 / 3)
        self.assertEqual(r["donor_support"], "mixed donor directions")

    def test_negative_mean_positive_majority(self):
        r = summarize_donors([-10, 1, 1])
        self.assertLess(r["mean_log2_contrast"], 0)
        self.assertEqual(r["agreement_with_mean"], 1 / 3)
        self.assertEqual(r["donor_support"], "mixed donor directions")

    def test_zeros_remain_in_denominator(self):
        r = summarize_donors([2, 2, 0])
        self.assertEqual(r["donor_support"], "2 of 3 positive")
        self.assertEqual(r["agreement_with_mean"], 2 / 3)
        self.assertEqual(summarize_donors([-2, 0])["donor_support"], "1 of 2 negative")

    def test_ties_and_zero_mean(self):
        r = summarize_donors([-1, -1, 1, 1])
        self.assertEqual(r["mean_direction"], "zero")
        self.assertEqual(r["donor_support"], "mixed donor directions")
        self.assertTrue(np.isnan(r["agreement_with_mean"]))

    def test_no_rounding_threshold_for_sign(self):
        self.assertEqual(summarize_donors([1e-12, 1e-12])["donor_support"], "2 of 2 positive")

    def test_nonfinite_rejected(self):
        for x in (np.nan, np.inf, -np.inf):
            with self.assertRaisesRegex(ValueError, "finite"):
                summarize_donors([1, x])

    def fixture(self):
        return pd.DataFrame([dict(scope="identity", cell_type="A", gene="g", normalization=m,
                                  replicate=d, log2_contrast=x)
                             for m, values in (("hvg_cpm", [-2, -1]), ("source_cpm", [2, -1]), ("tmm_cpm", [3, -1]))
                             for d, x in zip(("d1", "d2"), values)])

    def test_normalization_independent_of_donor_mixing(self):
        result = summarize_pairs(self.fixture())
        self.assertEqual(result.query("normalization == 'source_cpm'").donor_support.item(), "mixed donor directions")
        self.assertFalse(compare_methods(result).direction_stable.item())

    def test_duplicate_donor_rejected(self):
        pairs = self.fixture()
        with self.assertRaisesRegex(ValueError, "unique donor"):
            summarize_pairs(pd.concat([pairs, pairs.iloc[:1]]))

    def test_same_count_different_donors_rejected(self):
        pairs = self.fixture()
        pairs.loc[0, "replicate"] = "different"
        with self.assertRaisesRegex(ValueError, "same eligible"):
            compare_methods(summarize_pairs(pairs))

    def test_explicit_zero_pair_stratum(self):
        coverage = pd.DataFrame([dict(scope="identity", cell_type=x) for x in ("A", "absent")])
        result = summarize_pairs(self.fixture(), coverage)
        absent = result[result.cell_type.eq("absent")]
        self.assertTrue(absent.donor_support.eq("not estimable").all())
        compared = compare_methods(result).query("cell_type == 'absent'").iloc[0]
        self.assertTrue(pd.isna(compared.direction_stable))
        self.assertEqual(compared.normalization_status, "not estimable")

    def test_missing_condition_not_imputed(self):
        pb = pd.DataFrame([dict(replicate="d1", condition="ctrl", cell_type="A", n_cells=20, g__cpm=1)])
        pairs, summary = paired_directionality(pb, ["g"], ["replicate", "cell_type"], {})
        self.assertTrue(pairs.empty)
        self.assertEqual(summary.donor_support.item(), "not estimable")
        self.assertEqual(summary.n_pairs.item(), 0)

    def test_entirely_excluded_type_is_reported(self):
        pb = pd.DataFrame([dict(replicate="d1", condition=c, cell_type="A", n_cells=20, g__cpm=1) for c in ("ctrl", "stim")])
        _, summary = paired_directionality(pb, ["g"], ["replicate", "cell_type"], {}, pd.DataFrame({"cell_type": ["A", "absent"]}))
        self.assertEqual(summary.query("cell_type == 'absent'").donor_support.item(), "not estimable")


if __name__ == "__main__":
    unittest.main()
