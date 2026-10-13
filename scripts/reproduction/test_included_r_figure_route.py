#!/usr/bin/env python3
"""Safety and input-selection checks for the saved-results figure command."""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_selected_k6_figures as figures


class IncludedRFigureRouteTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.results = self.root / "data/processed/selected_model_k6/r"
        self.results.mkdir(parents=True)
        self.names = (
            "cogaps_K6_seed2_iter2000.metrics.json", "pattern_summary.csv",
            "pattern_activity_by_celltype_condition.csv",
            "pattern_activity_by_replicate_condition.csv",
            "pattern_cell_activities_with_metadata.csv.gz",
        )
        for name in self.names:
            (self.results / name).write_text("fixture\n")
        self.output = self.root / "learner_outputs/figures"
        self.bundle = self.root / "data/processed/selected_model_k6/directionality_source_wide"
        self.original = {name: getattr(figures, name) for name in (
            "ROOT", "RESULTS_DIR", "DIRECTION_GLOBAL", "SENSITIVITY_DIR", "IMG_DIR",
            "SOURCE_TABLE_DIR", "LEGACY_IMG_DIR", "METRICS_JSON", "PATTERN_SUMMARY",
            "CELLTYPE_CONDITION", "REPLICATE_CONDITION", "CELL_ACTIVITIES",
        )}
        figures.ROOT = self.root
        self.bundle_check = patch.object(figures, "check_bundle").start()
        self.contrast_check = patch.object(figures, "load_canonical_directionality").start()

    def tearDown(self):
        patch.stopall()
        for name, value in self.original.items():
            setattr(figures, name, value)
        self.temp.cleanup()

    def test_separate_inputs_and_new_output_only(self):
        figures.configure_included_r_route(self.output)
        self.bundle_check.assert_called_once_with(self.bundle)
        self.contrast_check.assert_called_once_with()
        self.assertEqual(figures.PATTERN_SUMMARY.parent, self.results)
        self.assertEqual(figures.DIRECTION_GLOBAL.parent, self.bundle / "r")
        self.assertEqual(figures.SENSITIVITY_DIR, self.bundle / "sensitivity")
        self.assertEqual(figures.IMG_DIR, self.output)
        self.assertEqual(figures.SOURCE_TABLE_DIR, self.output / "source_tables")
        self.assertIsNone(figures.LEGACY_IMG_DIR)
        self.assertTrue(self.output.is_dir())

    def test_existing_output_is_preserved(self):
        self.output.mkdir(parents=True)
        marker = self.output / "keep.txt"
        marker.write_text("unchanged")
        with self.assertRaisesRegex(ValueError, "already exists"):
            figures.configure_included_r_route(self.output)
        self.assertEqual(marker.read_text(), "unchanged")
        self.bundle_check.assert_not_called()

    def test_existing_empty_directory_is_preserved(self):
        self.output.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, "already exists"):
            figures.configure_included_r_route(self.output)

    def test_existing_file_is_preserved(self):
        target = self.root / "keep.txt"
        target.write_text("unchanged")
        with self.assertRaisesRegex(ValueError, "already exists"):
            figures.configure_included_r_route(target)
        self.assertEqual(target.read_text(), "unchanged")

    def test_broken_symlink_is_rejected(self):
        target = self.root / "broken"
        target.symlink_to(self.root / "absent")
        with self.assertRaisesRegex(ValueError, "already exists"):
            figures.configure_included_r_route(target)
        self.assertTrue(target.is_symlink())

    def test_project_inputs_and_published_outputs_are_protected(self):
        for name in ("data", "img", "docs", "scripts", "resources", "_extensions", ".git"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "outside"):
                figures.configure_included_r_route(self.root / name / "new_figures")

    def test_symlink_into_input_folder_is_rejected(self):
        alias = self.root / "input_alias"
        alias.symlink_to(self.root / "data", target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "outside"):
            figures.configure_included_r_route(alias / "new_figures")

    def test_missing_and_empty_inputs_do_not_create_output(self):
        for name in self.names:
            for empty in (False, True):
                with self.subTest(name=name, empty=empty):
                    path = self.results / name
                    path.write_text("") if empty else path.unlink()
                    with self.assertRaisesRegex(ValueError, "Missing or empty"):
                        figures.configure_included_r_route(self.output)
                    self.assertFalse(self.output.exists())
                    path.write_text("fixture\n")

    def test_bundle_failure_does_not_create_output(self):
        self.bundle_check.side_effect = ValueError("SHA-256 mismatch")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            figures.configure_included_r_route(self.output)
        self.contrast_check.assert_not_called()
        self.assertFalse(self.output.exists())

    def test_contrast_failure_does_not_create_output(self):
        self.contrast_check.side_effect = ValueError("Stimulation values disagree")
        with self.assertRaisesRegex(ValueError, "disagree"):
            figures.configure_included_r_route(self.output)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
