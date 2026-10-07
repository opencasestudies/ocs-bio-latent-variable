"""Test wrapper guards with fake Docker only; never start a container or model."""

import json
import os
from pathlib import Path
import pty
import shutil
import subprocess
import sys
import tempfile
import unittest


WRAPPER = Path(__file__).with_name("run_selected_model_r.sh")


class RWrapperTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="cs4 wrapper test ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve() / "project with spaces"
        self.wrapper = self.root / "scripts/reproduction/run_selected_model_r.sh"
        self.wrapper.parent.mkdir(parents=True)
        shutil.copy2(WRAPPER, self.wrapper)
        self.outdir = self.root / "data/processed/selected_model_k6/r"
        self.result = self.outdir / "cogaps_K6_seed2_iter2000.rds"
        self.bin = Path(self.temporary.name) / "bin"
        self.bin.mkdir()
        self.calls_path = Path(self.temporary.name) / "docker_calls.jsonl"
        self.bash = shutil.which("bash")
        self.assertIsNotNone(self.bash)
        # A restricted PATH prevents fallback to any real Docker executable.
        for name in ("dirname", "mkdir"):
            tool = shutil.which(name)
            self.assertIsNotNone(tool)
            (self.bin / name).symlink_to(tool)
        self.env = dict(os.environ)
        for name in ("FORCE_RERUN", "COGAPS_RUNTIME_IMAGE", "BASH_ENV", "ENV"):
            self.env.pop(name, None)
        self.env.update(PATH=str(self.bin), TEST_DOCKER_LOG=str(self.calls_path))

    def fake_docker(self):
        executable = self.bin / "docker"
        executable.write_text(
            f"#!{sys.executable}\n"
            "import json, os, sys\n"
            "with open(os.environ['TEST_DOCKER_LOG'], 'a') as log:\n"
            "    log.write(json.dumps(sys.argv[1:]) + '\\n')\n"
            "if sys.argv[1:] == ['info']:\n"
            "    if os.environ.get('TEST_INFO_EXIT', '0') != '0':\n"
            "        print('FAKE daemon connection failure', file=sys.stderr)\n"
            "    sys.exit(int(os.environ.get('TEST_INFO_EXIT', '0')))\n"
            "if sys.argv[1] == 'run':\n"
            "    print('FAKE Docker: model command recorded, not executed')\n"
            "    sys.exit(int(os.environ.get('TEST_RUN_EXIT', '0')))\n"
            "sys.exit(99)\n"
        )
        executable.chmod(0o755)

    def existing_result(self, content=b"existing model sentinel"):
        self.outdir.mkdir(parents=True)
        self.result.write_bytes(content)
        return content

    def run_wrapper(self, **environment):
        return subprocess.run(
            [self.bash, str(self.wrapper)],
            cwd=self.temporary.name,
            env={**self.env, **environment},
            capture_output=True,
            text=True,
            timeout=15,
        )

    def calls(self):
        if not self.calls_path.exists():
            return []
        return [json.loads(line) for line in self.calls_path.read_text().splitlines()]

    def test_missing_docker_explains_host_terminal_without_writing_outputs(self):
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("Docker is not available", result.stderr)
        self.assertIn("not in the RStudio Terminal", result.stderr)
        self.assertIn("check that Docker is installed", result.stderr)
        self.assertEqual(self.calls(), [])
        self.assertFalse(self.outdir.exists())

    def test_unreachable_daemon_preserves_error_and_never_calls_run(self):
        self.fake_docker()
        result = self.run_wrapper(TEST_INFO_EXIT="37")
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("FAKE daemon connection failure", result.stderr)
        self.assertIn("no model was started", result.stderr)
        self.assertIn("retry docker info", result.stderr)
        self.assertEqual(self.calls(), [["info"]])
        self.assertFalse(self.outdir.exists())

    def test_existing_result_blocks_before_docker(self):
        self.fake_docker()
        original = self.existing_result()
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("FORCE_RERUN=1", result.stderr)
        self.assertEqual(self.calls(), [])
        self.assertEqual(self.result.read_bytes(), original)

    def test_existing_empty_result_is_still_protected(self):
        self.existing_result(b"")
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertEqual(self.result.read_bytes(), b"")
        self.assertEqual(self.calls(), [])

    def test_existing_result_guard_precedes_missing_docker(self):
        original = self.existing_result()
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("Full R model already exists", result.stderr)
        self.assertEqual(self.result.read_bytes(), original)

    def test_only_exact_force_value_one_bypasses_overwrite_guard(self):
        self.fake_docker()
        original = self.existing_result()
        for value in ("0", "true", "yes", "2", ""):
            with self.subTest(force=value):
                result = self.run_wrapper(FORCE_RERUN=value)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(self.calls(), [])
                self.assertEqual(self.result.read_bytes(), original)

    def test_force_does_not_bypass_daemon_guard(self):
        self.fake_docker()
        original = self.existing_result()
        result = self.run_wrapper(FORCE_RERUN="1", TEST_INFO_EXIT="1")
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(self.calls(), [["info"]])
        self.assertEqual(self.result.read_bytes(), original)

    def test_ready_path_preserves_mount_image_and_model_settings(self):
        self.fake_docker()
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual(calls[0], ["info"])
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[1][:-1], [
            "run", "--rm", "--platform", "linux/amd64",
            "-v", f"{self.root}:/workspace/case-study",
            "-w", "/workspace/case-study",
            "othomas2/pycogaps-runtime-guide:0.3.1", "bash", "-lc",
        ])
        for setting in (
            "OMP_NUM_THREADS=4", "CoGAPS::compiledWithOpenMPSupport()",
            "Rscript scripts/cogaps_run_one_singleprocess_r.R",
            "--k 6", "--seed 2", "--n-iter 2000", "--use-sparse-opt",
            "--n-snapshots 10", "--take-pump-samples", "--force-rerun",
        ):
            self.assertIn(setting, calls[1][-1])
        self.assertTrue(self.outdir.is_dir())
        self.assertFalse(self.result.exists())

    def test_explicit_force_and_image_override_reach_fake_docker_only(self):
        self.fake_docker()
        original = self.existing_result()
        result = self.run_wrapper(FORCE_RERUN="1", COGAPS_RUNTIME_IMAGE="test-image:local")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls()[0], ["info"])
        self.assertIn("test-image:local", self.calls()[1])
        self.assertEqual(self.result.read_bytes(), original)

    def test_interactive_terminal_keeps_docker_tty_option(self):
        self.fake_docker()
        master, slave = pty.openpty()
        try:
            result = subprocess.run(
                [self.bash, str(self.wrapper)], cwd=self.temporary.name,
                env=self.env, stdout=slave, stderr=subprocess.PIPE,
                text=True, timeout=15,
            )
        finally:
            os.close(slave)
            os.close(master)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls()[1][:5], ["run", "--rm", "-t", "--platform", "linux/amd64"])
        self.assertFalse(self.result.exists())

    def test_container_failure_exit_status_is_preserved(self):
        self.fake_docker()
        result = self.run_wrapper(TEST_RUN_EXIT="42")
        self.assertEqual(result.returncode, 42, result.stderr)
        self.assertEqual(len(self.calls()), 2)
        self.assertFalse(self.result.exists())


if __name__ == "__main__":
    unittest.main()
