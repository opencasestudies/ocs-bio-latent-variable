"""Focused, dependency-free tests for the route-aware inventory."""

import hashlib
import io
import json
from contextlib import redirect_stdout
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from check_reproduction_inputs import ROUTES, check_route, format_summary, inspect_file, main


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

    def test_summary_groups_missing_input_and_existing_outputs(self):
        self.put("prepared")
        self.put("config")
        text = format_summary(self.check("preprocess"), self.root)
        self.assertEqual(text, "\n".join([
            "Stage: Rebuild preprocessing",
            "Paths below are relative to your project folder unless shown as absolute paths.",
            "", "Missing required input:", "  data/source",
            "", "Files this stage will replace:", "  data/prepared", "  data/config",
            "", "Not ready. Add the missing input and run this check again.",
            "This check has not changed any files.",
        ]))
        self.assertNotIn("data/dense", text)
        self.assertNotIn("python_model", text)

    def test_ready_summary_hides_unneeded_files_and_new_outputs(self):
        self.put("prepared")
        text = format_summary(self.check("fit-r"), self.root)
        self.assertEqual(text, "Stage: Fit the R model\n\n"
                         "Ready: the required files passed the checks.\n"
                         "This check has not changed any files.")

    def test_invalid_file_reasons_remain_visible(self):
        for content, reason in ((b"", "file is empty"),
                                (b"short", "byte size mismatch"),
                                (b"x" * self.specs["source"]["size_bytes"], "SHA-256 mismatch")):
            with self.subTest(reason=reason):
                self.put("source", content=content)
                text = format_summary(self.check("preprocess"), self.root)
                self.assertIn("Invalid required inputs:\n  data/source", text)
                self.assertIn(reason, text)
                self.assertIn("Not ready.", text)

    def test_unreadable_input_reason_remains_visible(self):
        self.put("source")
        with patch.object(Path, "open", side_effect=PermissionError("denied")):
            report = self.check("preprocess")
        self.assertIn("cannot read file: denied", format_summary(report, self.root))

    def test_verified_archive_copy_has_concise_placement_guidance(self):
        archived = self.put("source", archive=True)
        text = format_summary(self.check("preprocess"), self.root)
        self.assertIn(f"Verified download: {archived}", text)
        self.assertIn("Copy it to the project path above", text)
        self.assertIn("Not ready.", text)

    def test_invalid_archive_copy_is_still_reported(self):
        self.put("source", archive=True, content=b"bad")
        text = format_summary(self.check("preprocess"), self.root)
        self.assertIn("Download check failed:", text)
        self.assertIn("byte size mismatch", text)
        self.assertNotIn("Verified download:", text)

    def test_summary_preserves_model_collision_for_both_languages(self):
        for key in ("prepared", "dense", "r_model", "python_model"):
            self.put(key)
        for route, model in (("fit-r", "r_model"), ("fit-python", "python_model")):
            with self.subTest(route=route):
                text = format_summary(self.check(route), self.root)
                self.assertIn(f"Existing model files block fitting:\n  data/{model}", text)
                self.assertIn("refuses to overwrite", text)
                self.assertNotIn("Files this stage will replace:", text)
                self.assertIn("Not ready.", text)

    def test_replacement_warning_does_not_change_ready_result(self):
        for key in ("source", "prepared", "config"):
            self.put(key)
        report = self.check("preprocess")
        text = format_summary(report, self.root)
        self.assertTrue(report["ready"])
        self.assertIn("Files this stage will replace:", text)
        self.assertIn("Ready: the required files passed the checks.", text)

    def test_formatting_is_read_only_and_does_not_recheck_files(self):
        self.put("prepared")
        report = self.check("preprocess")
        before = deepcopy(report)
        with patch.object(Path, "exists", side_effect=AssertionError("unexpected recheck")), \
             patch.object(Path, "open", side_effect=AssertionError("unexpected file access")):
            format_summary(report, self.root)
        self.assertEqual(report, before)

    def test_cli_json_and_exit_codes_match_reports_for_all_routes(self):
        self.root = self.root.resolve()
        self.archive = self.archive.resolve()
        self.root.mkdir(exist_ok=True)
        (self.root / "index.qmd").touch()
        (self.root / "scripts").mkdir()
        (self.root / "data").mkdir()
        for populated in (False, True):
            if populated:
                for key in self.specs:
                    self.put(key)
            for route in ROUTES:
                with self.subTest(populated=populated, route=route):
                    report = self.check(route)
                    before = {str(p): p.read_bytes() for p in self.root.rglob("*") if p.is_file()}
                    args = ["--route", route, "--project-root", str(self.root), "--archive", str(self.archive)]
                    with patch("check_reproduction_inputs.read_specs", return_value=self.specs):
                        for as_json in (False, True):
                            output = io.StringIO()
                            with redirect_stdout(output):
                                code = main(args + (["--json"] if as_json else []))
                            self.assertEqual(code, 0 if report["ready"] else 1)
                            if as_json:
                                self.assertEqual(json.loads(output.getvalue()), report)
                            else:
                                self.assertEqual(output.getvalue(), format_summary(report, self.root) + "\n")
                    after = {str(p): p.read_bytes() for p in self.root.rglob("*") if p.is_file()}
                    self.assertEqual(after, before)


if __name__ == "__main__":
    unittest.main()
