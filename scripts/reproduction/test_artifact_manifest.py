"""Small inventory fixtures test failures without full models or scientific packages."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from check_artifact_manifest import ARCHIVE_INDEX, BUNDLE_RELATIVE, FIGURE_MANIFEST, check_manifest


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.manifest_path = self.root / "data/artifact_manifest.json"
        self.write(BUNDLE_RELATIVE / "README.md", "Included results\n")
        self.write(BUNDLE_RELATIVE / "r/result.csv", "value\n1\n")
        self.write(BUNDLE_RELATIVE / "bundle_manifest.json", json.dumps({"files": [
            {**self.record(BUNDLE_RELATIVE / "r/result.csv"), "path": "r/result.csv"}]}))
        self.write("scripts/builder.py", "# Current figure builder\n")
        self.write(FIGURE_MANIFEST, json.dumps({"validated_generator": self.record("scripts/builder.py"),
                                              "validated_against": [], "outputs": []}))
        self.write(ARCHIVE_INDEX, json.dumps({"record_url": "https://example.org/record/1", "doi": "example",
            "files": [{"target_path": "data/source/source.h5ad", "upload_filename": "source.h5ad",
                       "size_bytes": 8, "sha256": "a" * 64}]}))
        included = [BUNDLE_RELATIVE / name for name in ("README.md", "bundle_manifest.json", "r/result.csv")]
        included.append(FIGURE_MANIFEST)
        self.manifest = {"schema_version": "1.0", "included_files": [
            {**self.record(p), "github_required": True} for p in included],
            "case_study_model": {"directionality_result_dir": BUNDLE_RELATIVE.as_posix()},
            "external_archive_index": self.record(ARCHIVE_INDEX),
            "optional_external_or_local_files": [{"recommended_path": "data/source/source.h5ad",
                "upload_filename": "source.h5ad", "size_bytes": 8, "sha256": "a" * 64,
                "record_url": "https://example.org/record/1", "doi": "example", "github_required": False}]}

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def record(self, name):
        path = self.root / name
        return {"path": str(name), "size_bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    def check(self):
        self.manifest_path.write_text(json.dumps(self.manifest))
        return check_manifest(self.root)

    def test_optional_download_not_required(self):
        self.assertTrue(self.check()["ready"])
        self.assertFalse((self.root / "data/source/source.h5ad").exists())

    def test_missing_included(self):
        (self.root / BUNDLE_RELATIVE / "README.md").unlink()
        with self.assertRaisesRegex(ValueError, "Missing file"):
            self.check()

    def test_same_size_corruption(self):
        self.write(BUNDLE_RELATIVE / "r/result.csv", "value\n2\n")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            self.check()

    def test_empty(self):
        self.write(BUNDLE_RELATIVE / "README.md", "")
        self.manifest["included_files"][0].update(self.record(BUNDLE_RELATIVE / "README.md"))
        with self.assertRaisesRegex(ValueError, "Invalid included"):
            self.check()

    def test_duplicate(self):
        self.manifest["included_files"].append(self.manifest["included_files"][0])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.check()

    def test_path_escape(self):
        for path in ("../outside", "/tmp/outside"):
            with self.subTest(path=path):
                self.manifest["included_files"][0]["path"] = path
                with self.assertRaisesRegex(ValueError, "relative path"):
                    self.check()

    def test_symlink_escape(self):
        path = self.root / BUNDLE_RELATIVE / "README.md"
        path.unlink()
        path.symlink_to(self.root.parent / "outside")
        with self.assertRaisesRegex(ValueError, "leaves the bundle"):
            self.check()

    def test_staged(self):
        self.manifest["staged_artifact_bundles"] = [{"path": "old"}]
        with self.assertRaisesRegex(ValueError, "staged"):
            self.check()

    def test_missing_bundle_entry(self):
        self.manifest["included_files"].pop(2)
        with self.assertRaisesRegex(ValueError, "omits a directionality"):
            self.check()

    def test_optional_identity_mismatch(self):
        self.manifest["optional_external_or_local_files"][0]["upload_filename"] = "wrong.h5ad"
        with self.assertRaisesRegex(ValueError, "disagrees"):
            self.check()

    def test_duplicate_optional(self):
        self.manifest["optional_external_or_local_files"] *= 2
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.check()

    def test_local_presence_flag(self):
        self.manifest["optional_external_or_local_files"][0]["present_locally"] = True
        with self.assertRaisesRegex(ValueError, "author"):
            self.check()

    def test_stale_figure_builder(self):
        self.write("scripts/builder.py", "# A different builder\n")
        with self.assertRaisesRegex(ValueError, "mismatch"):
            self.check()

    def test_retired_file(self):
        self.write("data/results/pattern_direction_summary.csv", "old\n")
        with self.assertRaisesRegex(ValueError, "Superseded"):
            self.check()


if __name__ == "__main__":
    unittest.main()
