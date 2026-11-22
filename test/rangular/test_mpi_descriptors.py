"""Exercise root-owned descriptor input through real MPI broadcasts and MCP output."""
from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

SUITE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SUITE / "test"))
sys.path.insert(0, str(SUITE / "scripts"))
from rdensity.regression import fixture, floats, records
from test_csf_descriptors import Record, Shell, basic_blocks, write_fixture
import mcp_compare

ARGUMENTS: argparse.Namespace


class MpiDescriptorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="rgmpi-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        parquet_path, header_path, self.text_path = write_fixture(
            self.root / "inputs", ["1s"], ["2s", "2p-", "2p", "3s", "3p-", "3p"], basic_blocks()
        )
        self.binary_path = parquet_path

    def execute(
        self, label: str, ranks: int, options: list[str], *,
        interactive_input: str = "n\n", default_text: bool = False,
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        directory = self.root / label
        directory.mkdir()
        if default_text:
            shutil.copyfile(self.binary_path, directory / "rcsf_descriptors.parquet")
            shutil.copyfile(self.binary_path.parent / "csfs_header.toml", directory / "rcsf_header.toml")
        return self.run_mpi(ARGUMENTS.rangular, directory, ranks, options, interactive_input), directory

    @staticmethod
    def run_mpi(
        executable: Path, directory: Path, ranks: int, options: list[str], stdin: str
    ) -> subprocess.CompletedProcess[str]:
        environment = os.environ | {
            "MPI_TMP": str(directory), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "OMPI_ALLOW_RUN_AS_ROOT": "1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
            "OMPI_MCA_rmaps_base_oversubscribe": "1",
        }
        return subprocess.run(
            [str(ARGUMENTS.mpiexec), ARGUMENTS.numproc_flag, str(ranks), str(executable), *options],
            cwd=directory, input=stdin, text=True, capture_output=True,
            env=environment, timeout=30, check=False,
        )

    def successful_run(
        self, label: str, ranks: int, options: list[str], *,
        interactive_input: str = "n\n", default_text: bool = False,
    ) -> Path:
        result, directory = self.execute(
            label, ranks, options, interactive_input=interactive_input, default_text=default_text
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Execution complete", result.stdout)
        return directory

    @staticmethod
    def outputs(directory: Path) -> dict[str, bytes]:
        return {str(path.relative_to(directory)): path.read_bytes() for path in directory.glob("[0-9][0-9][0-9]/mcp.*")}

    def dataset(self, directory: Path) -> mcp_compare.Dataset:
        return mcp_compare.Dataset(directory)

    def test_full_interaction_and_empty_ranks(self) -> None:
        reference_dataset = None
        for ranks in (1, 2, 4):
            with self.subTest(ranks=ranks):
                text_directory = self.successful_run(f"text-{ranks}", ranks, [], default_text=True)
                binary_directory = self.successful_run(
                    f"binary-{ranks}", ranks,
                    ["--csf-descriptors", str(self.binary_path), "--timings", "timings.csv"],
                )
                self.assertFalse((binary_directory / "rcsf.inp").exists())
                self.assertEqual(len(self.outputs(text_directory)), ranks)
                input_comparison = mcp_compare.compare(self.dataset(text_directory), self.dataset(binary_directory),
                                                     atol=0.0, rtol=0.0)
                self.assertEqual(input_comparison["status"], "equal", input_comparison)
                for rank_index in range(ranks):
                    profile = binary_directory / f"{rank_index:03d}" / f"timings.rank{rank_index:03d}.csv"
                    with profile.open() as stream:
                        rows = list(csv.DictReader(stream))
                    self.assertTrue(all(float(row["wall_seconds"]) >= 0 for row in rows))
                    loads = [row for row in rows if row["stage"] == "csf_load_and_broadcast"]
                    self.assertEqual([int(row["records"]) for row in loads], [5, 2, 3])
                    self.assertTrue(all(row["storage"] == "descriptors" for row in loads))
                    self.assertEqual(rows[-1]["stage"], "program_before_finalize_including_interaction_input")
                current_dataset = self.dataset(binary_directory)
                if reference_dataset is None:
                    reference_dataset = self.dataset(text_directory)
                comparison = mcp_compare.compare(reference_dataset, current_dataset, atol=0.0, rtol=0.0)
                self.assertEqual(comparison["status"], "equal", comparison)
                self.assertGreater(sum(channel["compared"] for channel in comparison["channels"]), 0)

    def test_first_order_and_absolute_profile_paths(self) -> None:
        parquet_path, _, _ = write_fixture(
            self.root / "first-inputs", ["1s"], ["2s", "2p-", "2p", "3s", "3p-", "3p"],
            basic_blocks(), zero_order_counts=[1, 2, 2],
        )
        full = self.successful_run(
            "full", 2, ["--csf-descriptors", str(parquet_path)], interactive_input="n\n"
        )
        candidate = self.successful_run(
            "first", 2,
            ["--csf-descriptors", str(parquet_path), "--timings", str(self.root / "profile.csv")],
            interactive_input="y\n",
        )
        self.assertNotEqual(self.outputs(full), self.outputs(candidate))
        summary = (candidate / "000" / "fort.24").read_text()
        self.assertEqual(
            [line.strip() for line in summary.splitlines() if "zero-order space" in line],
            [f"CSFs 1--{count} constitute the zero-order space." for count in [1, 2, 2]],
        )
        self.assertTrue((self.root / "profile.rank000.csv").is_file())
        self.assertTrue((self.root / "profile.rank001.csv").is_file())

    def test_invalid_or_missing_zero_first_counts_stop_mpi_job(self) -> None:
        for case, counts, message in (
            ("missing", None, "missing block_info.zero_order_counts"),
            ("empty-prefix", [0, 1, 1], "positive zero_order_counts"),
            ("wrong-length", [1], "length does not match block count"),
            ("outside-block", [1, 3, 1], "invalid zero-order CSF count for block 2"),
        ):
            with self.subTest(case=case):
                parquet, _, _ = write_fixture(
                    self.root / f"inputs-{case}", ["1s"], ["2s", "2p-", "2p", "3s", "3p-", "3p"],
                    basic_blocks(), zero_order_counts=counts,
                )
                result, directory = self.execute(
                    case, 2, ["--csf-descriptors", str(parquet)], interactive_input="y\n"
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertFalse(self.outputs(directory))

    def test_rci_zero_first_prompt_is_independent_of_default_settings(self) -> None:
        if ARGUMENTS.rci is None:
            self.skipTest("rci_mpi is needed for the shared partition test")
        energies: dict[str, float] = {}
        for name, prefix, answer, counts in [
            ("default-full", "y\ncase\n", "n", None),
            ("default-full-with-counts", "y\ncase\n", "n", [1]),
            ("default-partition", "y\ncase\n", "y", [1]),
            ("nondefault-partition", "n\ncase\nn\nn\n", "y", [1]),
        ]:
            with self.subTest(case=name):
                directory = self.root / name
                fixture(directory, True, True)
                # Two first-order CSFs with nonzero Coulomb coupling. The
                # partition removes their off-diagonal term from the CI matrix.
                parquet, header, _ = write_fixture(
                    directory, [], ["1s", "3s", "2p"],
                    [[
                        Record((Shell("1s", 2),), 0),
                        Record((Shell("3s", 2),), 0),
                        Record((Shell("2p", 2, 0, 0),), 0),
                    ]],
                    zero_order_counts=counts,
                )
                parquet.rename(directory / "case_descriptors.parquet")
                header.rename(directory / "case_header.toml")
                for filename in ("case.c", "case.cm", "rcsf.inp"):
                    (directory / filename).unlink()
                result = self.run_mpi(
                    ARGUMENTS.rci, directory, 2, [],
                    prefix + answer + "\nn\nn\nn\nn\nn\nn\n1\n",
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertTrue((directory / "case.cimix").is_file())
                summary = (directory / "case.csum").read_text()
                self.assertEqual("CSFs 1--1 constitute the zero-order space" in summary, answer == "y")
                shift, shifted_energy = floats(records(directory / "case.cimix")[-2])
                energy = shift + shifted_energy
                self.assertTrue(math.isfinite(energy), f"{name}: nonfinite CI energy {energy}")
                energies[name] = energy
        self.assertEqual(len(energies), 4, f"RCI cases with finite energies: {list(energies)}")
        self.assertAlmostEqual(energies["default-full"], energies["default-full-with-counts"])
        self.assertNotAlmostEqual(energies["default-full-with-counts"], energies["default-partition"])
        self.assertAlmostEqual(energies["default-partition"], energies["nondefault-partition"])

    def test_relative_inputs_use_serial_directory(self) -> None:
        # Relative paths are resolved before access from the rank-specific cwd.
        relative_binary = "../inputs/csfs_descriptors.parquet"
        reference = self.successful_run("relative-text", 2, ["--csf-descriptors", str(self.binary_path)])
        candidate = self.successful_run("relative-binary", 2, ["--csf-descriptors", relative_binary])
        comparison = mcp_compare.compare(self.dataset(reference), self.dataset(candidate), atol=0.0, rtol=0.0)
        self.assertEqual(comparison["status"], "equal", comparison)

    def test_invalid_descriptor_stops_mpi_job(self) -> None:
        import pyarrow as arrow
        from pyarrow import parquet
        for label in ("truncated", "bad-state"):
            with self.subTest(case=label):
                path = self.root / f"{label}_descriptors.parquet"
                shutil.copyfile(self.binary_path.parent / "csfs_header.toml", self.root / f"{label}_header.toml")
                if label == "truncated":
                    path.write_bytes(self.binary_path.read_bytes()[:-1])
                else:
                    table = parquet.read_table(self.binary_path)
                    columns = [c.to_pylist() for c in table.columns]
                    columns[1][0] = 1  # invalid full-shell seniority
                    parquet.write_table(arrow.Table.from_arrays(
                        [arrow.array(v, type=arrow.int8()) for v in columns], schema=table.schema), path)
                result, _ = self.execute(label, 2, ["--csf-descriptors", str(path)])
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("CSF descriptor", result.stderr)

    def test_existing_profile_is_preserved(self) -> None:
        profile_path = self.root / "existing.rank000.csv"
        profile_path.write_text("keep this file\n")
        result, _ = self.execute(
            "existing-profile", 2,
            ["--csf-descriptors", str(self.binary_path), "--timings", str(self.root / "existing.csv")],
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(profile_path.read_text(), "keep this file\n")

    def test_command_line_validation_and_help(self) -> None:
        for options in (["--csf-descriptors"], ["--unknown"],
                        ["--csf-text", "input.c", "--csf-descriptors", "input.bin"],
                        ["--csf-descriptors", "mcp000.30"]):
            result = subprocess.run([str(ARGUMENTS.rangular), *options], cwd=self.root,
                                    text=True, capture_output=True, timeout=10, check=False)
            self.assertNotEqual(result.returncode, 0)
        result = subprocess.run([str(ARGUMENTS.rangular), "--help"], cwd=self.root,
                                text=True, capture_output=True, timeout=10, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--csf-descriptors", result.stdout)
        self.assertIn("rank000.csv", result.stdout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rangular", type=lambda value: Path(value).resolve(), required=True)
    parser.add_argument("--rci", type=lambda value: Path(value).resolve())
    parser.add_argument("--mpiexec", type=lambda value: Path(value).resolve(), required=True)
    parser.add_argument("--numproc-flag", default="-np")
    ARGUMENTS, remaining = parser.parse_known_args()
    unittest.main(argv=[sys.argv[0], *remaining])
