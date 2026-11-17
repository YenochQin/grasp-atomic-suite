"""Exercise root-owned descriptor input through real MPI broadcasts and MCP output."""
from __future__ import annotations

import argparse
import csv
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
from test_csf_descriptors import basic_blocks, write_fixture
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
        interactive_input: str = "y\n", default_text: bool = False,
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        directory = self.root / label
        directory.mkdir()
        if default_text:
            shutil.copyfile(self.binary_path, directory / "rcsf_descriptors.parquet")
            shutil.copyfile(self.binary_path.parent / "csfs_header.toml", directory / "rcsf_header.toml")
        environment = os.environ | {
            "MPI_TMP": str(directory), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "OMPI_ALLOW_RUN_AS_ROOT": "1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
            "OMPI_MCA_rmaps_base_oversubscribe": "1",
        }
        result = subprocess.run(
            [str(ARGUMENTS.mpiexec), ARGUMENTS.numproc_flag, str(ranks), str(ARGUMENTS.rangular), *options],
            cwd=directory, input=interactive_input, text=True, capture_output=True,
            env=environment, timeout=30, check=False,
        )
        return result, directory

    def successful_run(
        self, label: str, ranks: int, options: list[str], *,
        interactive_input: str = "y\n", default_text: bool = False,
    ) -> Path:
        result, directory = self.execute(
            label, ranks, options, interactive_input=interactive_input, default_text=default_text
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Execution complete", result.stdout)
        return directory

    @staticmethod
    def outputs(directory: Path) -> dict[str, bytes]:
        return {path.name: path.read_bytes() for path in directory.glob("[0-9][0-9][0-9]/mcp*.*")}

    def dataset(self, directory: Path) -> mcp_compare.Dataset:
        view = directory / "view"
        view.mkdir(exist_ok=True)
        shutil.copyfile(self.text_path, view / "rcsf.inp")
        for path in directory.glob("[0-9][0-9][0-9]/mcp*.*"):
            (view / path.name).symlink_to(path)
        return mcp_compare.Dataset(view)

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
                self.assertEqual(len(self.outputs(text_directory)), 6 * ranks)
                self.assertEqual(self.outputs(text_directory), self.outputs(binary_directory))
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
        # MPI GETINF requires 1 < ICCUT < block length, unlike the serial prompt.
        parquet_path, header_path, text_path = write_fixture(
            self.root / "first-inputs", ["1s"], ["2s", "2p-", "2p", "3s", "3p-", "3p"],
            [basic_blocks()[0]],
        )
        binary_path = parquet_path
        reference = self.successful_run(
            "first-text", 2, ["--csf-descriptors", str(parquet_path)], interactive_input="n\n2\n"
        )
        candidate = self.successful_run(
            "first-binary", 2,
            ["--csf-descriptors", str(binary_path), "--timings", str(self.root / "profile.csv")],
            interactive_input="n\n2\n",
        )
        self.assertEqual(self.outputs(reference), self.outputs(candidate))
        self.assertTrue((self.root / "profile.rank000.csv").is_file())
        self.assertTrue((self.root / "profile.rank001.csv").is_file())

    def test_relative_inputs_use_serial_directory(self) -> None:
        # Relative paths are resolved before access from the rank-specific cwd.
        relative_binary = "../inputs/csfs_descriptors.parquet"
        reference = self.successful_run("relative-text", 2, ["--csf-descriptors", str(self.binary_path)])
        candidate = self.successful_run("relative-binary", 2, ["--csf-descriptors", relative_binary])
        self.assertEqual(self.outputs(reference), self.outputs(candidate))

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
                    columns[2][0] = 1  # invalid full-shell seniority: reaches deep Fortran validation
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
    parser.add_argument("--mpiexec", type=lambda value: Path(value).resolve(), required=True)
    parser.add_argument("--numproc-flag", default="-np")
    ARGUMENTS, remaining = parser.parse_known_args()
    unittest.main(argv=[sys.argv[0], *remaining])
