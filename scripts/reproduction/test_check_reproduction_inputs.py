"""Focused, dependency-free tests for the route-aware inventory."""

import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from check_reproduction_inputs import check_route, inspect_file


class InventoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "project"
        self.root.mkdir()
        self.archive = Path(self.temp.name) / "archive"
        self.archive.mkdir()
        self.specs = {}
        for key in ("source", "prepared", "config", "dense", "r_model", "checkpoint", "python_model"):
            content = (key + " reference").encode()
            self.specs[key] = {"target_path": f"data/{key}", "upload_filename": f"{key}.download",
                               "size_bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
        for key in ("r_summary", "python_summary"):
            self.specs[key] = {"target_path": f"data/{key}"}

    def put(self, key, archive=False, content=None):
        spec = self.specs[key]
        path = self.archive / spec["upload_filename"] if archive else self.root / spec["target_path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((key + " reference").encode() if content is None else content)
        return path

    def check(self, route, generated=()):
        return check_route(self.root, self.archive, route, generated, self.specs)

    def row(self, report, key):
        return next(row for row in report["files"] if row["file"] == key)

    def test_empty_archive_is_not_ready_for_preprocessing(self):
        report = self.check("preprocess")
        self.assertFalse(report["ready"])
        self.assertEqual(report["archive"]["status"], "archive directory empty")
        self.assertEqual(self.row(report, "source")["status"], "missing required input")
        self.assertEqual(self.row(report, "prepared")["status"], "output to be generated")

    def test_learner_route_does_not_require_archive_or_full_models(self):
        for key in ("prepared", "config", "r_summary"):
            self.put(key)
        report = self.check("learner-r")
        self.assertTrue(report["ready"])
        self.assertEqual(self.row(report, "python_model")["status"], "not needed for this route")

    def test_fitting_r_needs_prepared_data_not_saved_results(self):
        self.put("prepared")
        report = self.check("fit-r")
        self.assertTrue(report["ready"])
        self.assertEqual(self.row(report, "r_model")["status"], "output to be generated")
        self.assertEqual(self.row(report, "source")["status"], "not needed for this route")

    def test_export_r_requires_full_r_not_python(self):
        self.put("prepared")
        self.assertFalse(self.check("export-r")["ready"])
        self.put("r_model")
        self.assertTrue(self.check("export-r")["ready"])

    def test_python_fit_needs_dense_input_not_saved_python_model(self):
        self.put("prepared")
        self.assertFalse(self.check("fit-python")["ready"])
        self.put("dense")
        self.assertTrue(self.check("fit-python")["ready"])

    def test_python_inspection_requires_its_full_model(self):
        self.put("prepared")
        self.put("python_model")
        self.assertTrue(self.check("inspect-python")["ready"])

    def test_zero_byte_file_is_invalid(self):
        self.put("source", content=b"")
        row = self.row(self.check("preprocess"), "source")
        self.assertEqual(row["status"], "invalid required input")
        self.assertIn("empty", row["detail"])

    def test_same_size_corruption_fails_sha256(self):
        self.put("source", content=b"x" * self.specs["source"]["size_bytes"])
        self.assertIn("SHA-256 mismatch", self.row(self.check("preprocess"), "source")["detail"])

    def test_wrong_size_is_invalid(self):
        self.put("source", content=b"wrong size")
        self.assertIn("byte size mismatch", self.row(self.check("preprocess"), "source")["detail"])

    def test_directory_is_not_an_input_file(self):
        path = self.root / self.specs["source"]["target_path"]
        path.mkdir(parents=True)
        self.assertIn("not a regular file", self.row(self.check("preprocess"), "source")["detail"])

    def test_unreadable_file_is_invalid(self):
        path = self.put("source")
        with patch.object(Path, "open", side_effect=PermissionError("denied")):
            ok, detail = inspect_file(path, self.specs["source"])
        self.assertFalse(ok)
        self.assertIn("cannot read", detail)

    def test_archive_copy_is_verified_but_not_automatically_used_or_copied(self):
        self.put("source", archive=True)
        before = list(self.root.rglob("*"))
        report = self.check("preprocess")
        self.assertFalse(report["ready"])
        self.assertIn("verified archive copy", self.row(report, "source")["detail"])
        self.assertEqual(before, list(self.root.rglob("*")))

    def test_corrupt_archive_copy_is_not_advertised_as_verified(self):
        self.put("source", archive=True, content=b"bad")
        row = self.row(self.check("preprocess"), "source")
        self.assertNotIn("verified archive copy", row["detail"])
        self.assertIn("byte size mismatch", row["archive_check"])

    def test_declared_generated_input_is_not_compared_to_original_model(self):
        self.put("prepared")
        self.put("r_model", content=b"a different validated generated model")
        self.assertFalse(self.check("export-r")["ready"])
        report = self.check("export-r", ["r_model"])
        self.assertTrue(report["ready"])
        self.assertIn("declared regenerated", self.row(report, "r_model")["detail"])

    def test_generated_missing_or_empty_still_fails(self):
        self.assertFalse(self.check("fit-r", ["prepared"])["ready"])
        self.put("prepared", content=b"")
        self.assertFalse(self.check("fit-r", ["prepared"])["ready"])

    def test_source_cannot_be_declared_generated(self):
        with self.assertRaisesRegex(ValueError, "never exempt"):
            self.check("preprocess", ["source"])

    def test_unused_generated_flag_is_rejected(self):
        with self.assertRaises(ValueError):
            self.check("fit-r", ["python_model"])

    def test_existing_full_result_stops_fitting_without_overwriting(self):
        self.put("prepared")
        path = self.put("r_model")
        before = path.read_bytes()
        report = self.check("fit-r")
        self.assertFalse(report["ready"])
        self.assertIn("refuses to overwrite", self.row(report, "r_model")["detail"])
        self.assertEqual(before, path.read_bytes())


if __name__ == "__main__":
    unittest.main()
