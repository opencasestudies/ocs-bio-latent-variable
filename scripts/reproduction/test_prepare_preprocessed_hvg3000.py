"""Run with python -m unittest discover -s scripts/reproduction -p 'test_*.py'."""

import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from prepare_preprocessed_hvg3000 import build_cogaps_input, ensure_condition_column


class ConditionTests(unittest.TestCase):
    def data(self, **columns):
        return ad.AnnData(
            sparse.csr_matrix([[1, 2, 0], [0, 3, 4]], dtype=np.float32),
            obs=pd.DataFrame(columns, index=["cell1", "cell2"]),
            var=pd.DataFrame(index=["g1", "g2", "g3"]),
        )

    def test_copies_label_without_changing_source_or_order(self):
        data = self.data(label=pd.Categorical(["ctrl", "stim"]), replicate=["d1", "d1"])
        before = data.obs.copy(deep=True)
        self.assertEqual(ensure_condition_column(data), "label")
        pd.testing.assert_frame_equal(data.obs[before.columns], before)
        pd.testing.assert_series_equal(data.obs.condition, before.label, check_names=False)

    def test_existing_condition_is_authoritative(self):
        data = self.data(condition=["ctrl", "stim"], label=["stim", "ctrl"])
        before = data.obs.copy(deep=True)
        self.assertEqual(ensure_condition_column(data), "condition")
        pd.testing.assert_frame_equal(data.obs, before)

    def test_missing_columns_fail(self):
        with self.assertRaisesRegex(ValueError, "Missing obs"):
            ensure_condition_column(self.data())

    def test_missing_values_fail_without_fallback(self):
        with self.assertRaisesRegex(ValueError, "missing treatment"):
            ensure_condition_column(self.data(condition=["ctrl", None], label=["ctrl", "stim"]))

    def test_unknown_labels_fail(self):
        with self.assertRaisesRegex(ValueError, "no other labels"):
            ensure_condition_column(self.data(label=["control", "stim"]))

    def test_single_condition_fails(self):
        with self.assertRaisesRegex(ValueError, "both"):
            ensure_condition_column(self.data(label=["ctrl", "ctrl"]))

    def test_dense_orientation_preserves_condition_and_values(self):
        data = self.data(label=["ctrl", "stim"])
        ensure_condition_column(data)
        result = build_cogaps_input(data)
        np.testing.assert_array_equal(result.X, data.X.toarray().T)
        self.assertEqual(result.X.dtype, np.float64)
        pd.testing.assert_frame_equal(result.var, data.obs)
        self.assertEqual(list(result.obs_names), list(data.var_names))


if __name__ == "__main__":
    unittest.main()
