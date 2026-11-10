"""Exercise the runner with real loaders and controlled failing/different producers."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from test_csf_descriptors import basic_blocks, write_fixture
from test_mcp_compare import TERMS
from test_mcp_compare import fixture as mcp_fixture

SUITE_ROOT = Path(__file__).resolve().parents[1]
RUNNER = SUITE_ROOT / "scripts/compare_rangular_versions.py"
ARGUMENTS: argparse.Namespace


class ComparisonRunnerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="rangular comparison ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.parquet_path, self.header_path, self.csf_path = write_fixture(
            self.root / "source inputs",
            ["1s"],
            ["2s", "2p-", "2p", "3s", "3p-", "3p"],
            basic_blocks(),
        )
        self.output = self.root / "results"

    def run_runner(
        self,
        *additional: str,
        reference: Path | None = None,
        candidate: Path | None = None,
        csf_path: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(RUNNER),
                "--reference",
                str(reference or ARGUMENTS.rangular),
                "--candidate",
                str(candidate or ARGUMENTS.rangular),
                "--csf",
                str(csf_path or self.csf_path),
                "--output",
                str(self.output),
                *additional,
            ],
            text=True,
            capture_output=True,
            timeout=45,
            check=False,
            env=os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
        )

    def report(self) -> dict:
        return json.loads((self.output / "summary.json").read_text())

    def executable(self, name: str, body: str) -> Path:
        path = self.root / name
        path.write_text("#!/usr/bin/env python3\n" + body)
        path.chmod(0o755)
        return path

    def producer(self, name: str, data: Path) -> Path:
        return self.executable(
            name,
            f"""from pathlib import Path
import shutil
import sys
for source in Path({str(data)!r}).glob('mcp.*'):
    shutil.copyfile(source, Path.cwd() / source.name)
if '--timings' in sys.argv:
    Path(sys.argv[sys.argv.index('--timings')+1]).write_text(
        'stage,block,channel,wall_seconds,records,storage\\n'
        'program_total_including_interaction_input,0,-1,0.001,-1,unspecified\\n')
""",
        )

    def test_real_three_paths_and_repeated_execution(self) -> None:
        result = self.run_runner(
            "--descriptors",
            str(self.parquet_path),
            "--header",
            str(self.header_path),
            "--verifier",
            str(ARGUMENTS.verifier),
            "--repeat",
            "2",
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = self.report()
        self.assertEqual(report["status"], "equal")
        self.assertTrue(report["same_executable_bytes"])
        self.assertEqual(len(report["runs"]), 6)
        self.assertEqual(len(report["comparisons"]), 5)
        self.assertEqual(report["array_verification"]["status"], "verified")
        self.assertEqual(report["preparation"]["rows"], 10)
        self.assertTrue(
            all(
                comparison.get("byte_identical", True)
                for comparison in report["comparisons"]
            )
        )
        self.assertTrue(all(run["wall_seconds"] > 0 for run in report["runs"]))
        self.assertEqual(report["runs"][3]["mode"], "candidate_text")
        for run in report["runs"]:
            self.assertTrue((Path(run["directory"]) / "stdout.log").exists())
            if run["mode"] == "candidate_descriptors":
                self.assertTrue(run["text_provenance_added_after_execution"])
        previous_report = (self.output / "summary.json").read_bytes()
        repeated = self.run_runner()
        self.assertEqual(repeated.returncode, 2)
        self.assertEqual(previous_report, (self.output / "summary.json").read_bytes())

    def test_prepared_binary_and_first_order_transcript(self) -> None:
        from prepare_csf_descriptors import prepare_descriptors

        binary_path = self.root / "prepared input.bin"
        prepare_descriptors(self.parquet_path, self.header_path, binary_path)
        stdin_path = self.root / "answers.txt"
        stdin_path.write_text("n\n1\n1\n1\n")
        original_digest = hashlib.sha256(self.csf_path.read_bytes()).hexdigest()
        result = self.run_runner(
            "--binary",
            str(binary_path),
            "--verifier",
            str(ARGUMENTS.verifier),
            "--stdin-file",
            str(stdin_path),
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = self.report()
        self.assertNotIn("preparation", report)
        self.assertEqual(
            original_digest, hashlib.sha256(self.csf_path.read_bytes()).hexdigest()
        )
        self.assertTrue(
            all(run["mcp_inspection"]["header"]["lfordr"] for run in report["runs"])
        )

    def test_coefficient_difference_is_failure_and_has_no_speed_ratio(self) -> None:
        reference_data = self.root / "reference_data"
        candidate_data = self.root / "candidate_data"
        mcp_fixture(reference_data)
        altered_terms = [(*TERMS[0][:-1], 0.75), *TERMS[1:]]
        mcp_fixture(candidate_data, terms=altered_terms)
        result = self.run_runner(
            reference=self.producer("reference", reference_data),
            candidate=self.producer("candidate", candidate_data),
            csf_path=reference_data / "rcsf.inp",
        )
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        report = self.report()
        self.assertEqual(report["status"], "different")
        self.assertNotIn("reference_over_candidate_wall_ratio", report)
        comparison = json.loads(Path(report["comparisons"][0]["report"]).read_text())
        self.assertGreater(comparison["difference_count"], 0)

    def test_zero_exit_without_mcp_is_not_success(self) -> None:
        result = self.run_runner(
            reference=self.executable("early_stop", "raise SystemExit(0)\n")
        )
        self.assertEqual(result.returncode, 2)
        self.assertEqual(self.report()["status"], "error")
        self.assertIn("no mcp", self.report()["error"])
        self.assertTrue((self.output / "trial_001_reference_text/stdout.log").exists())

    def test_timeout_preserves_failure_report(self) -> None:
        result = self.run_runner(
            "--timeout",
            "0.1",
            reference=self.executable("sleeping", "import time\ntime.sleep(30)\n"),
        )
        self.assertEqual(result.returncode, 2)
        self.assertEqual(self.report()["error_type"], "TimeoutExpired")
        execution = json.loads(
            (self.output / "trial_001_reference_text/execution.json").read_text()
        )
        self.assertEqual(execution["status"], "timeout")

    def test_verifier_legacy_stop_does_not_start_calculations(self) -> None:
        fake_binary = self.root / "input.bin"
        fake_binary.write_bytes(b"not used by deliberately broken verifier")
        verifier = self.executable("broken_verifier", "raise SystemExit(0)\n")
        result = self.run_runner(
            "--binary", str(fake_binary), "--verifier", str(verifier)
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("without its success marker", self.report()["error"])
        self.assertFalse(list(self.output.glob("trial_*")))

    def test_invalid_options_fail_before_output_creation(self) -> None:
        for options in (
            ("--descriptors", str(self.parquet_path)),
            ("--repeat", "0"),
            ("--timeout", "nan"),
        ):
            with self.subTest(options=options):
                result = self.run_runner(*options)
                self.assertEqual(result.returncode, 2)
                self.assertFalse(self.output.exists())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rangular", type=lambda value: Path(value).resolve(), required=True
    )
    parser.add_argument(
        "--verifier", type=lambda value: Path(value).resolve(), required=True
    )
    ARGUMENTS, remaining = parser.parse_known_args()
    unittest.main(argv=[sys.argv[0], *remaining])
