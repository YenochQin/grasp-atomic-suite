"""Exact serial/OpenMP MCP regression with actual concurrent angular calls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from prepare_csf_descriptors import prepare_descriptors
from test_csf_descriptors import Record, Shell, basic_blocks, write_fixture

ARGUMENTS: argparse.Namespace


class AngularOpenMPTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="rangular omp ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def run_angular(
        self,
        name: str,
        text_path: Path,
        *,
        threads: int | None = None,
        binary_path: Path | None = None,
        chunk_size: int = 256,
        memory_mib: int = 64,
        answers: str = "y\n",
    ) -> tuple[dict[str, str], str, list[dict[str, str]]]:
        directory = self.root / name
        directory.mkdir()
        executable = ARGUMENTS.serial if threads is None else ARGUMENTS.parallel
        command = [str(executable), "--timings", "timings.csv"]
        if binary_path is None:
            shutil.copyfile(text_path, directory / "rcsf.inp")
        else:
            command += ["--csf-descriptors", str(binary_path)]
        if threads is not None:
            command += [
                "--threads",
                str(threads),
                "--pair-chunk",
                str(chunk_size),
                "--pair-buffer-mib",
                str(memory_mib),
            ]
        result = subprocess.run(
            command,
            cwd=directory,
            input=answers,
            text=True,
            capture_output=True,
            timeout=90,
            check=False,
            env=os.environ
            | {
                "OMP_NUM_THREADS": "2",
                "OMP_DYNAMIC": "true",
                "OPENBLAS_NUM_THREADS": "1",
                "OMP_PROC_BIND": "false",
            },
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Execution complete", result.stdout)
        hashes = {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in directory.glob("mcp.*")
        }
        self.assertTrue(hashes)
        with (directory / "timings.csv").open() as stream:
            timings = list(csv.DictReader(stream))
        self.assertEqual(
            timings[-1]["stage"], "program_total_including_interaction_input"
        )
        return hashes, result.stdout, timings

    def test_concurrent_columns_and_block_reset(self) -> None:
        peel = [f"{principal_number}s" for principal_number in range(2, 22)]
        singlets = [Record((Shell(label, 2),), 0) for label in peel]
        pairs = [
            (left, right)
            for left_index, left in enumerate(peel)
            for right in peel[left_index + 1 :]
        ]
        singlets += [
            Record((Shell(left, 1, 1, 1), Shell(right, 1, 1, 1)), 0)
            for left, right in pairs
        ]
        triplets = [
            Record((Shell(left, 1, 1, 1), Shell(right, 1, 1, 1)), 2)
            for left, right in pairs
        ]
        parquet_path, header_path, text_path = write_fixture(
            self.root / "inputs", ["1s"], peel, [singlets, triplets]
        )
        binary_path = self.root / "input.bin"
        prepare_descriptors(parquet_path, header_path, binary_path)
        for interaction, answers in [("full", "y\n"), ("first", "n\n3\n3\n")]:
            expected, _, _ = self.run_angular(
                f"serial_{interaction}", text_path, answers=answers
            )
            for threads in (1, 2, 4, 8):
                with self.subTest(interaction=interaction, threads=threads):
                    actual, stdout, timings = self.run_angular(
                        f"omp_{interaction}_{threads}",
                        text_path,
                        threads=threads,
                        binary_path=binary_path if threads % 4 == 0 else None,
                        chunk_size=7,
                        answers=answers,
                    )
                    self.assertEqual(actual, expected)
                    batches = [
                        row
                        for row in timings
                        if row["stage"] == "parallel_batch_compute"
                    ]
                    self.assertEqual(len(batches), 2)
                    self.assertTrue(
                        all(row["storage"] == "omp-memory" for row in batches)
                    )
                    if threads > 1:
                        active = [
                            int(value)
                            for value in re.findall(r"active workers=(\d+)", stdout)
                        ]
                        self.assertTrue(any(value > 1 for value in active), stdout)

    def test_spill_and_memory_transfer_are_identical(self) -> None:
        _, _, text_path = write_fixture(
            self.root / "inputs",
            ["1s"],
            ["2s", "2p-", "2p", "3s", "3p-", "3p"],
            basic_blocks(),
        )
        expected, _, _ = self.run_angular("serial", text_path)
        for threads in (1, 4):
            actual, stdout, timings = self.run_angular(
                f"spill_{threads}",
                text_path,
                threads=threads,
                chunk_size=3,
                memory_mib=0,
            )
            self.assertEqual(actual, expected)
            self.assertTrue(any(row["storage"] == "omp-spill" for row in timings))
            self.assertTrue(
                any(
                    int(value) > 0
                    for value in re.findall(r"spilled chunks=(\d+)", stdout)
                )
            )

    def test_open_f_and_high_j_shared_tables(self) -> None:
        fixtures = Path(__file__).resolve().parents[1] / "jj2lsj"
        for fixture_name in ("even", "odd", "f3"):
            text_path = fixtures / f"{fixture_name}.c"
            expected, _, _ = self.run_angular(f"serial_{fixture_name}", text_path)
            for threads in (2, 8):
                with self.subTest(fixture=fixture_name, threads=threads):
                    actual, _, _ = self.run_angular(
                        f"omp_{fixture_name}_{threads}",
                        text_path,
                        threads=threads,
                        chunk_size=1,
                    )
                    self.assertEqual(actual, expected)

    def test_four_active_shells_and_seniority(self) -> None:
        records = []
        for outer_number in (3, 4):
            for intermediate_two_j in (0, 2):
                records.append(
                    Record(
                        (
                            Shell("2s", 1, 1, 1),
                            Shell("2p-", 1, 1, 1, intermediate_two_j),
                            Shell(f"{outer_number}s", 1, 1, 1, 1),
                            Shell(f"{outer_number}p-", 1, 1, 1),
                        ),
                        0,
                    )
                )
        records += [
            Record((Shell("2s", 2), Shell("2p-", 2)), 0),
            Record((Shell("3s", 2), Shell("3p-", 2)), 0),
        ]
        f_states = [
            Record((Shell("4f", 4, 4, seniority, printed_seniority=seniority),), 4)
            for seniority in (2, 4)
        ]
        for case_name, peel, blocks in [
            ("four_shells", ["2s", "2p-", "3s", "3p-", "4s", "4p-"], [records]),
            ("seniority", ["4f"], [f_states]),
        ]:
            _, _, text_path = write_fixture(self.root / case_name, [], peel, blocks)
            expected, _, _ = self.run_angular(f"serial_{case_name}", text_path)
            for threads in (1, 4, 8):
                actual, _, _ = self.run_angular(
                    f"omp_{case_name}_{threads}",
                    text_path,
                    threads=threads,
                    chunk_size=1,
                )
                self.assertEqual(actual, expected)

    def test_parallel_options_and_serial_guard(self) -> None:
        for executable, options in [
            (ARGUMENTS.serial, ["--threads", "2"]),
            (ARGUMENTS.parallel, ["--threads", "0"]),
            (ARGUMENTS.parallel, ["--pair-chunk", "0"]),
            (ARGUMENTS.parallel, ["--pair-buffer-mib", "-1"]),
        ]:
            result = subprocess.run(
                [str(executable), *options],
                cwd=self.root,
                input="",
                text=True,
                capture_output=True,
                check=False,
                timeout=10,
            )
            self.assertNotEqual(result.returncode, 0)

    def test_help_preserves_deferred_length_argument(self) -> None:
        for executable in (ARGUMENTS.serial, ARGUMENTS.parallel):
            result = subprocess.run(
                [str(executable), "--help"],
                text=True,
                capture_output=True,
                check=False,
                timeout=10,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("--csf-descriptors", result.stdout)

    def test_comparison_runner_sets_candidate_threads(self) -> None:
        _, _, text_path = write_fixture(
            self.root / "inputs",
            ["1s"],
            ["2s", "2p-", "2p", "3s", "3p-", "3p"],
            basic_blocks(),
        )
        runner = (
            Path(__file__).resolve().parents[2] / "scripts/compare_rangular_versions.py"
        )
        output = self.root / "comparison"
        result = subprocess.run(
            [
                sys.executable,
                str(runner),
                "--reference",
                str(ARGUMENTS.serial),
                "--candidate",
                str(ARGUMENTS.parallel),
                "--candidate-threads",
                "4",
                "--csf",
                str(text_path),
                "--output",
                str(output),
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        summary = json.loads((output / "summary.json").read_text())
        self.assertEqual(summary["candidate_threads"], 4)
        for run in summary["runs"]:
            expected = "1" if run["mode"] == "reference_text" else "4"
            self.assertEqual(run["thread_environment"]["OMP_NUM_THREADS"], expected)
            self.assertEqual(run["thread_environment"]["OPENBLAS_NUM_THREADS"], "1")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--serial", type=lambda value: Path(value).resolve(), required=True
    )
    parser.add_argument(
        "--parallel", type=lambda value: Path(value).resolve(), required=True
    )
    ARGUMENTS, remaining = parser.parse_known_args()
    unittest.main(argv=[sys.argv[0], *remaining])
