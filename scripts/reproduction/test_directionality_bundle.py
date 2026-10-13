"""Regression checks for standalone, saved-result imports."""

import csv
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from check_directionality_bundle import BUNDLE_RELATIVE, check_bundle, contained, validate_report


SOURCE = Path(__file__).resolve().parents[2] / BUNDLE_RELATIVE


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle = self.root / "bundle"
        shutil.copytree(SOURCE, self.bundle)
        self.report = self.bundle / "sensitivity/r_gene_report_global.csv"

    def edit_manifest(self, change):
        path = self.bundle / "bundle_manifest.json"
        manifest = json.loads(path.read_text())
        change(manifest)
        path.write_text(json.dumps(manifest))

    def edit_report(self, change):
        with self.report.open(newline="") as stream:
            reader = csv.DictReader(stream)
            fields, rows = reader.fieldnames, list(reader)
        change(fields, rows)
        with self.report.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def test_bundle_without_models_source_or_packages(self):
        result = check_bundle(self.bundle)
        self.assertTrue(result["ready"])
        self.assertFalse(result["prepared_identity_checked"])
        self.assertEqual(sum(result["canonical_report_rows"].values()), 3060)

    def test_missing_file(self):
        self.report.unlink()
        with self.assertRaisesRegex(ValueError, "Missing file"):
            check_bundle(self.bundle)

    def test_same_size_corruption(self):
        content = self.report.read_bytes()
        self.report.write_bytes(b"X" + content[1:])
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            check_bundle(self.bundle)

    def test_empty_file(self):
        self.report.write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "byte-size"):
            check_bundle(self.bundle)

    def test_old_schema(self):
        self.edit_manifest(lambda manifest: manifest.update(reporting_version="old"))
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            check_bundle(self.bundle)

    def test_unindexed_component(self):
        self.edit_manifest(lambda manifest: manifest.update(files=[r for r in manifest["files"] if r["path"] != "r/directionality_method.json"]))
        with self.assertRaisesRegex(ValueError, "Missing indexed"):
            check_bundle(self.bundle)

    def test_missing_required_report_field(self):
        def remove(fields, rows):
            fields.remove("n_zero")
            for row in rows:
                del row["n_zero"]
        self.edit_report(remove)
        with self.assertRaisesRegex(ValueError, "Missing donor"):
            validate_report(self.report)

    def test_confidence_column_rejected(self):
        def legacy(fields, rows):
            fields.append("evidence_strength")
            for row in rows:
                row["evidence_strength"] = "strong"
        self.edit_report(legacy)
        with self.assertRaisesRegex(ValueError, "Legacy confidence"):
            validate_report(self.report)

    def test_donor_counts(self):
        self.edit_report(lambda fields, rows: rows[0].update(n_positive="999"))
        with self.assertRaisesRegex(ValueError, "donor counts"):
            validate_report(self.report)

    def test_wrong_sensitivity_flag(self):
        self.edit_report(lambda fields, rows: rows[0].update(normalization_status="normalization-sensitive mean sign"))
        with self.assertRaisesRegex(ValueError, "normalization flag"):
            validate_report(self.report)

    def test_wrong_agreement(self):
        self.edit_report(lambda fields, rows: rows[0].update(agreement_with_mean="0.5"))
        with self.assertRaisesRegex(ValueError, "agreement with mean"):
            validate_report(self.report)

    def test_wrong_prepared_file(self):
        path = self.root / "prepared.h5ad"
        path.write_text("not the prepared input")
        with self.assertRaisesRegex(ValueError, "byte-size"):
            check_bundle(self.bundle, path)

    def test_paths_stay_in_bundle(self):
        for name in ("../escape", "/tmp/escape"):
            with self.assertRaises(ValueError):
                contained(self.bundle, name)

    def test_symlink_cannot_leave_bundle(self):
        (self.bundle / "outside").symlink_to(self.root)
        with self.assertRaisesRegex(ValueError, "leaves the bundle"):
            contained(self.bundle, "outside/input.csv")


if __name__ == "__main__":
    unittest.main()
