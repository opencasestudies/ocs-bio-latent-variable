#!/usr/bin/env python3

import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_selected_k6_figures as figures
from directionality_normalization import file_record
from directionality_reporting import summarize_donors


class FigureContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        figures.RESULTS_DIR = self.root
        figures.SENSITIVITY_DIR = self.root
        figures.DIRECTION_GLOBAL = self.root / "pattern_gene_directionality_global.csv"
        common = dict(cell_type="A", gene="g", normalization="source_cpm", n_donors=2,
                      normalization_status="same mean sign across methods")
        pd.DataFrame([{**common, **summarize_donors([2, 2]), "scope": "identity"}]).to_csv(self.root / "identity_source_summary.csv", index=False)
        pd.DataFrame([{**common, **summarize_donors([0.2, 0.2]), "scope": "stimulation_global"}]).to_csv(self.root / "stimulation_source_summary.csv", index=False)
        pd.DataFrame([dict(pattern="Pattern1", gene="g", mean_log2fc_stim_vs_ctrl=0.2)]).to_csv(figures.DIRECTION_GLOBAL, index=False)
        self.normalization = {"source": {"sha256": "source"}, "prepared": {"sha256": "prepared"}}
        self.write_manifests()

    def tearDown(self):
        self.temp.cleanup()

    def write_manifests(self):
        method = {"schema_version": "directionality-sensitivity-v2", "reporting_version": "donor-support-v1", "pseudocount_CPM": 0.5,
                  "normalization": self.normalization,
                  "outputs": [file_record(self.root / name) for name in ("identity_source_summary.csv", "stimulation_source_summary.csv")]}
        (self.root / "sensitivity_method.json").write_text(json.dumps(method))
        (self.root / "directionality_method.json").write_text(json.dumps({"normalization": self.normalization,
            "outputs": [file_record(figures.DIRECTION_GLOBAL)]}))

    def test_saved_values_are_accepted(self):
        identity, direction = figures.load_canonical_directionality()
        self.assertEqual(identity.mean_log2_contrast.iloc[0], 2)
        self.assertEqual(direction.mean_log2fc_stim_vs_ctrl.iloc[0], 0.2)

    def test_no_missing_manifest_fallback(self):
        (self.root / "sensitivity_method.json").unlink()
        with self.assertRaisesRegex(ValueError, "fallback"):
            figures.load_canonical_directionality()

    def test_corrupted_table(self):
        with (self.root / "identity_source_summary.csv").open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(ValueError, "Checksum"):
            figures.load_canonical_directionality()

    def test_duplicate_identity(self):
        path = self.root / "identity_source_summary.csv"
        table = pd.read_csv(path)
        pd.concat([table, table]).to_csv(path, index=False)
        self.write_manifests()
        with self.assertRaisesRegex(ValueError, "duplicate"):
            figures.load_canonical_directionality()

    def test_missing_gene(self):
        path = self.root / "stimulation_source_summary.csv"
        table = pd.read_csv(path)
        table["gene"] = "other"
        table.to_csv(path, index=False)
        self.write_manifests()
        with self.assertRaisesRegex(ValueError, "disagree"):
            figures.load_canonical_directionality()

    def test_independent_stimulation_disagreement(self):
        table = pd.read_csv(figures.DIRECTION_GLOBAL)
        table["mean_log2fc_stim_vs_ctrl"] = -0.2
        table.to_csv(figures.DIRECTION_GLOBAL, index=False)
        self.write_manifests()
        with self.assertRaisesRegex(ValueError, "disagree"):
            figures.load_canonical_directionality()

    def test_changed_pseudocount(self):
        path = self.root / "sensitivity_method.json"
        method = json.loads(path.read_text())
        method["pseudocount_CPM"] = 1
        path.write_text(json.dumps(method))
        with self.assertRaisesRegex(ValueError, "contract"):
            figures.load_canonical_directionality()

    def test_mixed_source_provenance(self):
        path = self.root / "directionality_method.json"
        method = json.loads(path.read_text())
        method["normalization"]["source"]["sha256"] = "different"
        path.write_text(json.dumps(method))
        with self.assertRaisesRegex(ValueError, "provenance"):
            figures.load_canonical_directionality()

    def test_legacy_confidence_rejected(self):
        table = pd.read_csv(figures.DIRECTION_GLOBAL)
        table["evidence_strength"] = "strong"
        table.to_csv(figures.DIRECTION_GLOBAL, index=False)
        self.write_manifests()
        with self.assertRaisesRegex(ValueError, "Legacy"):
            figures.load_canonical_directionality()

    def test_donor_counts_required(self):
        path = self.root / "identity_source_summary.csv"
        pd.read_csv(path).drop(columns="n_positive").to_csv(path, index=False)
        self.write_manifests()
        with self.assertRaisesRegex(ValueError, "donor support"):
            figures.load_canonical_directionality()

    def test_negative_donor_counts_rejected(self):
        path = self.root / "identity_source_summary.csv"
        table = pd.read_csv(path)
        table["n_positive"], table["n_negative"] = -1, 3
        table.to_csv(path, index=False)
        self.write_manifests()
        with self.assertRaisesRegex(ValueError, "Invalid donor counts"):
            figures.load_canonical_directionality()

    def test_sensitive_bars_have_hatching_and_counts(self):
        import matplotlib.pyplot as plt
        table = pd.DataFrame([{**summarize_donors([2, -1]), "normalization_status": "normalization-sensitive mean sign"}])
        fig, ax = plt.subplots()
        figures.draw_reported_contrast(ax, [0], table, "mean_log2_contrast", ["red"])
        self.assertEqual(ax.patches[0].get_hatch(), "///")
        self.assertEqual(ax.texts[0].get_text(), "1 / 1 / 0")
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
