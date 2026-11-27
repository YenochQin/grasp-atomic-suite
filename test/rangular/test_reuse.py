"""Exact MCP regression for angular reuse, including incomplete radial families."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import shutil
import sys
import tempfile
import unittest
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "test"))
from test_csf_descriptors import Record, Shell, write_fixture

ARGUMENTS: argparse.Namespace


def reuse_fixture(root: Path, radial_count: int = 6) -> tuple[Path, Path, Path]:
    # Interleaving symmetries forbids grouping/sorting orbitals by kappa. The
    # subset of occupations is intentionally not a complete Cartesian family.
    s = [f"{n}s" for n in range(2, radial_count + 2)]
    p = [f"{n}p" for n in range(2, radial_count + 2)]
    peel = [label for pair in zip(s, p, strict=True) for label in pair]
    order = {label: index for index, label in enumerate(peel)}
    records: list[Record] = []
    for labels in combinations(s, 4):
        for intermediate in (0, 2):
            records.append(Record(tuple(
                Shell(label, 1, 1, 1, two_k=coupling)
                for label, coupling in zip(labels, (None, intermediate, 1, None), strict=True)
            ), 0))
    for labels in combinations(s, 2):
        records.append(Record(tuple(Shell(label, 2) for label in labels), 0))
    for closed in s:
        for left, right in combinations([label for label in s if label != closed], 2):
            # The fixture writer puts total J at the last printed shell;
            # trailing closed shells use a different legacy text layout.
            if order[closed] > order[right]:
                continue
            shells = [Shell(closed, 2), Shell(left, 1, 1, 1), Shell(right, 1, 1, 1)]
            records.append(Record(tuple(sorted(shells, key=lambda shell: order[shell.label])), 0))
    for left, right in combinations(p, 2):
        for two_j, seniority in ((0, 0), (4, 2)):
            records.append(Record((Shell(left, 2, two_j, seniority), Shell(right, 2, two_j, seniority)), 0))
        records.append(Record((Shell(left, 3, 3, 1), Shell(right, 1, 3, 1)), 0))
    # A second block with the same J/parity tests cache lifetime and block order.
    return write_fixture(root, ["1s"], peel, [records, records[::7]])


class AngularReuseTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="angular reuse ")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def compare(self, csf: Path, *options: str, memory: int = 8) -> dict:
        if "--descriptors" not in options:
            from native_fixture import descriptors
            fixture = self.root / "input.c"
            shutil.copyfile(csf, fixture)
            descriptors(fixture)
            options = (*options, "--descriptors", str(self.root / "input_desc.parquet"),
                       "--header", str(self.root / "input_header.toml"), "--verifier", str(ARGUMENTS.verifier))
        output = self.root / "comparison"
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/compare_rangular_versions.py"),
             "--reference", str(ARGUMENTS.rangular), "--candidate", str(ARGUMENTS.rangular),
             "--csf", str(csf), "--output", str(output),
             "--candidate-reuse-cache-mib", str(memory), "--reference-timings", *options],
            capture_output=True, text=True, timeout=120,
            env=os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads((output / "summary.json").read_text())
        self.assertEqual(report["status"], "equal")
        self.assertTrue(all(item.get("byte_identical", True) for item in report["comparisons"]))
        return report

    def test_default_and_explicit_descriptor_reuse_with_overlap_and_coupling_changes(self) -> None:
        parquet, header, csf = reuse_fixture(self.root / "input")
        report = self.compare(csf, "--descriptors", str(parquet), "--header", str(header),
                              "--verifier", str(ARGUMENTS.verifier), "--repeat", "2")
        stats = [run["angular_reuse"] for run in report["runs"] if run["mode"] != "reference_default"]
        self.assertTrue(all(item["totals"]["hits"] > 0 for item in stats))
        self.assertTrue(all(item["totals"]["misses"] > 0 for item in stats))
        self.assertTrue(all(item == stats[0] for item in stats))
        self.assertEqual(len(stats[0]["blocks"]), 2)
        self.assertTrue(all(block["misses"] > 0 for block in stats[0]["blocks"]))
        if sys.platform in {"darwin", "linux"}:
            self.assertTrue(all(run["peak_rss_bytes"] > 0 for run in report["runs"]))

    def test_first_order_subset_preserves_original_structure(self) -> None:
        _, _, csf = reuse_fixture(self.root / "input")
        answers = self.root / "stdin.txt"
        answers.write_text("n\n5\n3\n")
        report = self.compare(csf, "--stdin-file", str(answers))
        self.assertTrue(all(run["mcp_inspection"]["header"]["lfordr"] for run in report["runs"]))

    def test_full_cache_falls_back_without_changing_output(self) -> None:
        report = self.compare(ROOT / "data/test.c", memory=1)
        stats = next(run["angular_reuse"] for run in report["runs"] if run["mode"] == "candidate_default")
        self.assertGreater(stats["totals"]["not_stored"], 0)
        self.assertLessEqual(stats["peak_cache_reserved_bytes"], 1048576)

    def test_distinct_subshell_terms_are_not_merged(self) -> None:
        self.compare(ROOT / "test/jj2lsj/f3.c")

    def test_reject_invalid_or_unsupported_options_before_computation(self) -> None:
        for options in (
            ["--reuse-cache-mib", "-1"], ["--reuse-cache-mib", "4097"],
            ["--reuse-stats", "stats.csv"],
            ["--reuse-cache-mib", "1", "--reuse-stats", "case.mcp.30"],
            ["--reuse-cache-mib", "1", "--timings", "same.csv", "--reuse-stats", "same.csv"],
        ):
            with self.subTest(options=options):
                result = subprocess.run([str(ARGUMENTS.rangular), *options], cwd=self.root,
                                        input="case\ny\n", capture_output=True, text=True, timeout=10)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(list(self.root.glob("case.mcp.*")))
        existing = self.root / "stats.csv"
        existing.write_text("preserve me\n")
        result = subprocess.run([str(ARGUMENTS.rangular), "--reuse-cache-mib", "1", "--reuse-stats", str(existing)],
                                cwd=self.root, capture_output=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(existing.read_text(), "preserve me\n")
        if ARGUMENTS.parallel:
            result = subprocess.run([str(ARGUMENTS.parallel), "--reuse-cache-mib", "1"], cwd=self.root,
                                    capture_output=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rangular", type=Path, required=True)
    parser.add_argument("--verifier", type=Path, required=True)
    parser.add_argument("--parallel", type=Path)
    ARGUMENTS, remaining = parser.parse_known_args()
    ARGUMENTS.rangular = ARGUMENTS.rangular.resolve()
    ARGUMENTS.verifier = ARGUMENTS.verifier.resolve()
    if ARGUMENTS.parallel:
        ARGUMENTS.parallel = ARGUMENTS.parallel.resolve()
    unittest.main(argv=[sys.argv[0], *remaining])
