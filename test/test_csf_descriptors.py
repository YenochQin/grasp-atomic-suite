"""Container rejection, independent text-loader oracle, and rangular regressions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import prepare_csf_descriptors as prepare
import pyarrow as arrow
from pyarrow import parquet
from verify_csf_descriptors import verify_descriptors

ARGUMENTS: argparse.Namespace


@dataclass(frozen=True)
class Shell:
    label: str
    occupation: int
    two_j: int | None = None
    seniority: int = 0
    two_k: int | None = None
    printed_seniority: int | None = None


@dataclass(frozen=True)
class Record:
    shells: tuple[Shell, ...]
    total_two_j: int
    parity: int = 1


def angular_text(two_j: int) -> str:
    return str(two_j // 2) if two_j % 2 == 0 else f"{two_j}/2"


def text_record(record: Record) -> str:
    occupation_line = ""
    state_line = ""
    coupling_line = list(" " * (9 * len(record.shells) + 2))
    for shell_position, shell in enumerate(record.shells):
        display_label = shell.label if shell.label.endswith("-") else shell.label + " "
        occupation_line += f"{display_label:>5}({shell.occupation:2})"
        state_field = list(" " * 9)
        if shell.two_j is not None:
            value = angular_text(shell.two_j)
            state_field[9 - len(value) :] = value
            if shell.printed_seniority is not None:
                state_field[3:5] = f"{shell.printed_seniority};"
        state_line += "".join(state_field)
        if shell.two_k is not None:
            value = angular_text(shell.two_k)
            field_end = (shell_position + 1) * 9 + 2
            coupling_line[field_end + 1 - len(value) : field_end + 1] = value
    value = angular_text(record.total_two_j)
    field_end = 9 * len(record.shells)
    coupling_line[field_end + 1 - len(value) : field_end + 1] = value
    coupling_line[-1] = "+" if record.parity == 1 else "-"
    return (
        occupation_line
        + "\n"
        + state_line.rstrip()
        + "\n"
        + "".join(coupling_line)
        + "\n"
    )


def write_fixture(
    root: Path, core: list[str], peel: list[str], blocks: list[list[Record]]
) -> tuple[Path, Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    header_lines = [
        "Core subshells:",
        "  ".join(core),
        "Peel subshells:",
        "  ".join(peel),
        "CSF(s):",
    ]
    text_path = root / "rcsf.inp"
    text_path.write_text(
        "\n".join(header_lines)
        + "\n"
        + " *\n".join(
            "".join(text_record(record) for record in block) for block in blocks
        ),
        encoding="ascii",
    )
    header_path = root / "csfs_header.toml"
    header_path.write_text(
        "[header_info]\nheader_lines = "
        + json.dumps(header_lines)
        + "\n[block_info]\nblock_lengths = "
        + json.dumps([len(block) for block in blocks])
        + "\n"
    )
    names = [
        f"sub{shell_index}_{channel}"
        for shell_index in range(len(peel))
        for channel in ("n", "2j", "v", "2k")
    ] + ["total_two_j", "parity"]
    metadata = {
        "descriptor_version": "2",
        "descriptor_dtype": "int8",
        "max_two_j": "126",
        "seniority_encoding": "physical",
        "channels_per_subshell": "4",
        "subshell_count": str(len(peel)),
        "missing_sentinel": "-1",
        "normalized": "false",
        "source_header_sha256": hashlib.sha256(header_path.read_bytes()).hexdigest(),
        "peel_subshells": json.dumps(peel),
        "feature_columns": json.dumps(names[:-2]),
        "global_columns": json.dumps(names[-2:]),
    }
    values = []
    for block in blocks:
        for record in block:
            row = [0, -1, 0, -1] * len(peel) + [record.total_two_j, record.parity]
            for shell in record.shells:
                offset = 4 * peel.index(shell.label)
                row[offset : offset + 4] = [
                    shell.occupation,
                    -1 if shell.two_j is None else shell.two_j,
                    shell.seniority,
                    -1 if shell.two_k is None else shell.two_k,
                ]
            values.append(row)
    schema = arrow.schema(
        [arrow.field(name, arrow.int8(), nullable=False) for name in names],
        metadata=metadata,
    )
    table = arrow.Table.from_arrays(
        [
            arrow.array(column, type=arrow.int8())
            for column in zip(*values, strict=True)
        ],
        schema=schema,
    )
    parquet_path = root / "descriptors.parquet"
    parquet.write_table(
        table, parquet_path, row_group_size=2, compression="zstd", use_dictionary=True
    )
    return parquet_path, header_path, text_path


def basic_blocks() -> list[list[Record]]:
    return [
        [
            Record((Shell("2s", 2),), 0),
            Record((Shell("2p-", 2),), 0),
            Record((Shell("2p", 2, 0, 0),), 0),
            Record((Shell("3s", 2),), 0),
            Record((Shell("2s", 1, 1, 1), Shell("3s", 1, 1, 1)), 0),
        ],
        [
            Record((Shell("2p", 2, 4, 2),), 4),
            Record((Shell("2p-", 1, 1, 1), Shell("2p", 1, 3, 1)), 4),
        ],
        [
            Record((Shell("2s", 1, 1, 1), Shell("2p-", 1, 1, 1)), 2, -1),
            Record((Shell("2s", 1, 1, 1), Shell("2p", 1, 3, 1)), 2, -1),
            Record((Shell("3s", 1, 1, 1), Shell("3p-", 1, 1, 1)), 2, -1),
        ],
    ]


class DescriptorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.parquet_path, self.header_path, self.text_path = write_fixture(
            self.root, ["1s"], ["2s", "2p-", "2p", "3s", "3p-", "3p"], basic_blocks()
        )
        self.binary_path = self.root / "csfs.bin"

    def convert(self, **options: int) -> dict[str, object]:
        return prepare.prepare_descriptors(
            self.parquet_path, self.header_path, self.binary_path, **options
        )

    def native(self, binary: Path | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                str(ARGUMENTS.verifier),
                str(binary or self.binary_path),
                str(self.text_path),
            ],
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
        )

    def test_container_and_batch_boundaries(self) -> None:
        report = self.convert(batch_rows=1)
        self.assertEqual(report["rows"], 10)
        self.assertEqual(report["blocks"], 3)
        encoded = self.binary_path.read_bytes()
        self.assertEqual(encoded[:8], b"GCSDESC1")
        self.assertEqual(struct.unpack_from("<I", encoded, 40)[0], 4)
        self.assertEqual(len(encoded), 108 + 8 * 7 + 16 * 3 + 26 * 10)
        second_path = self.root / "larger_batches.bin"
        prepare.prepare_descriptors(
            self.parquet_path, self.header_path, second_path, batch_rows=7
        )
        self.assertEqual(encoded, second_path.read_bytes())
        with self.assertRaises(FileExistsError):
            self.convert()
        self.assertEqual(encoded, self.binary_path.read_bytes())

    def test_reject_metadata_and_header_mismatch(self) -> None:
        original = parquet.read_table(self.parquet_path)
        for key, value in [
            (b"descriptor_version", b"1"),
            (b"seniority_encoding", b"printed"),
            (b"descriptor_dtype", b"int32"),
            (b"missing_sentinel", b"0"),
            (b"source_header_sha256", b"0" * 64),
            (b"normalized", b"true"),
        ]:
            with self.subTest(key=key):
                metadata = dict(original.schema.metadata or {})
                metadata[key] = value
                parquet.write_table(
                    original.replace_schema_metadata(metadata), self.parquet_path
                )
                with self.assertRaises(ValueError):
                    self.convert()
                self.assertFalse(self.binary_path.exists())
        parquet.write_table(original, self.parquet_path)
        self.header_path.write_text(self.header_path.read_text() + "# changed\n")
        with self.assertRaisesRegex(ValueError, "source_header_sha256"):
            self.convert()

    def test_reject_schema_and_channels_without_leaving_output(self) -> None:
        original = parquet.read_table(self.parquet_path)
        for column_name, replacement in [
            ("sub0_n", -1),
            ("sub0_n", 3),
            ("sub0_v", -1),
            ("sub0_2j", 127),
            ("parity", 0),
            ("total_two_j", -1),
        ]:
            with self.subTest(column=column_name, replacement=replacement):
                column_index = original.schema.get_field_index(column_name)
                values = original.column(column_index).to_pylist()
                values[-1] = replacement
                changed = original.set_column(
                    column_index,
                    original.schema.field(column_index),
                    arrow.array(values, type=arrow.int8()),
                )
                parquet.write_table(changed, self.parquet_path)
                with self.assertRaises(ValueError):
                    self.convert(batch_rows=1)
                self.assertFalse(self.binary_path.exists())
        changed = original.cast(
            arrow.schema(
                [
                    arrow.field(name, arrow.int32(), nullable=False)
                    for name in original.column_names
                ],
                metadata=original.schema.metadata,
            )
        )
        parquet.write_table(changed, self.parquet_path)
        with self.assertRaisesRegex(ValueError, "Int8"):
            self.convert()

    def test_text_array_equality(self) -> None:
        self.convert()
        result = verify_descriptors(
            ARGUMENTS.verifier, self.binary_path, self.text_path
        )
        self.assertIn("Verified array values: 6350", result.stdout)

    def test_parquet_footer_metadata_without_arrow_schema_metadata(self) -> None:
        original = parquet.read_table(self.parquet_path)
        unbound_schema_table = original.replace_schema_metadata(None)
        with parquet.ParquetWriter(
            self.parquet_path, unbound_schema_table.schema
        ) as writer:
            writer.add_key_value_metadata(original.schema.metadata)
            writer.write_table(unbound_schema_table)
        self.assertIsNone(parquet.ParquetFile(self.parquet_path).schema_arrow.metadata)
        self.convert()
        verify_descriptors(ARGUMENTS.verifier, self.binary_path, self.text_path)

    def test_array_mismatch_is_localized(self) -> None:
        self.parquet_path, self.header_path, self.text_path = write_fixture(
            self.root,
            [],
            ["4f"],
            [[Record((Shell("4f", 4, 4, 2, printed_seniority=2),), 4)]],
        )
        self.convert()
        encoded = bytearray(self.binary_path.read_bytes())
        encoded[108 + 8 + 16 + 2] = 4
        self.binary_path.write_bytes(encoded)
        result = self.native()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("array mismatch", result.stderr)
        self.assertIn(
            "block, CSF, shell/slot, channel, text, descriptor: 1 1 1 1 2 4",
            result.stderr,
        )

    def test_upper_angular_byte_boundary(self) -> None:
        shells = [
            Shell(
                f"{principal_number}h",
                2,
                20,
                2,
                20 * (principal_number - 5) if principal_number > 6 else None,
            )
            for principal_number in range(6, 11)
        ]
        shells += [Shell("11g", 2, 16, 2, 116), Shell("12f", 2, 12, 2)]
        self.parquet_path, self.header_path, self.text_path = write_fixture(
            self.root,
            [],
            [shell.label for shell in shells],
            [[Record(tuple(shells), 126)]],
        )
        self.convert()
        verify_descriptors(ARGUMENTS.verifier, self.binary_path, self.text_path)

    def test_missing_single_shell_j_is_inferred(self) -> None:
        self.parquet_path, self.header_path, self.text_path = write_fixture(
            self.root, ["1s"], ["2p"], [[Record((Shell("2p", 1, 3, 1),), 3, -1)]]
        )
        self.convert()
        encoded = bytearray(self.binary_path.read_bytes())
        encoded[108 + 8 * 2 + 16 + 1] = 255
        self.binary_path.write_bytes(encoded)
        verify_descriptors(ARGUMENTS.verifier, self.binary_path, self.text_path)

    def test_open_zero_j_closed_shells_and_printed_zero_coupling(self) -> None:
        blocks = [
            [
                Record(
                    (Shell("2s", 1, 1, 1), Shell("2p", 2, 0, 0), Shell("3s", 1, 1, 1)),
                    0,
                ),
                Record(
                    (
                        Shell("2s", 1, 1, 1),
                        Shell("2p-", 1, 1, 1, 0),
                        Shell("3s", 1, 1, 1, 1),
                        Shell("3p-", 1, 1, 1),
                    ),
                    0,
                ),
                Record(
                    (Shell("2s", 1, 1, 1), Shell("2p-", 2), Shell("3s", 1, 1, 1)), 0
                ),
                Record(
                    (Shell("2s", 2), Shell("2p-", 1, 1, 1), Shell("3p-", 1, 1, 1)), 0
                ),
                Record(
                    (Shell("2s", 1, 1, 1), Shell("3s", 1, 1, 1, 0), Shell("3p-", 2)), 0
                ),
            ]
        ]
        self.parquet_path, self.header_path, self.text_path = write_fixture(
            self.root, ["1s"], ["2s", "2p-", "2p", "3s", "3p-"], blocks
        )
        self.convert()
        verify_descriptors(ARGUMENTS.verifier, self.binary_path, self.text_path)

    def test_seniority_and_adjacent_equal_symmetry_blocks(self) -> None:
        blocks = [
            [
                Record((Shell("4f", 4, 4, 2, printed_seniority=2),), 4),
                Record((Shell("4f", 4, 4, 4, printed_seniority=4),), 4),
            ],
            [Record((Shell("4f", 4, 4, 2, printed_seniority=2),), 4)],
        ]
        self.parquet_path, self.header_path, self.text_path = write_fixture(
            self.root, [], ["4f"], blocks
        )
        self.convert()
        verify_descriptors(ARGUMENTS.verifier, self.binary_path, self.text_path)

    def test_semantic_corruption_and_unsupported_states(self) -> None:
        self.convert()
        original = self.binary_path.read_bytes()
        payload_offset = 108 + 8 * 7 + 16 * 3
        mutations = [
            (payload_offset + 2, 1),  # closed shell seniority
            (payload_offset + 4 + 1, 0),  # empty shell printed J
            (payload_offset + 26 * 5 + 4 * 2 + 1, 2),
        ]  # impossible p^2 J=1
        for byte_offset, byte_value in mutations:
            with self.subTest(offset=byte_offset):
                encoded = bytearray(original)
                encoded[byte_offset] = byte_value
                self.binary_path.write_bytes(encoded)
                self.assertNotEqual(self.native().returncode, 0)
        # j=9/2, q=4, v=4, J=4 needs a w label absent from descriptors.
        self.parquet_path, self.header_path, self.text_path = write_fixture(
            self.root, [], ["5g"], [[Record((Shell("5g", 4, 8, 4),), 8)]]
        )
        self.binary_path.unlink()
        self.convert()
        result = subprocess.run(
            [str(ARGUMENTS.rangular), "--csf-descriptors", str(self.binary_path)],
            cwd=self.root,
            input="y\n",
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
            env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("additional w label", result.stderr)

    def test_binary_header_rejections(self) -> None:
        self.convert()
        original = self.binary_path.read_bytes()
        mutations = [original[:-1], original + b"x", b"WRONG!!!" + original[8:]]
        for byte_offset, value in [
            (8, 2),
            (12, 0x04030201),
            (16, 200),
            (24, 51),
            (28, 1),
            (32, 0),
        ]:
            encoded = bytearray(original)
            struct.pack_into("<I", encoded, byte_offset, value)
            mutations.append(bytes(encoded))
        for sequence, encoded in enumerate(mutations):
            with self.subTest(sequence=sequence):
                self.binary_path.write_bytes(encoded)
                self.assertNotEqual(self.native().returncode, 0)

    def test_verify_wrapper_rejects_legacy_stop_without_success(self) -> None:
        self.convert()
        self.text_path.write_text("not a CSF file\n")
        with self.assertRaisesRegex(ValueError, "verification failed"):
            verify_descriptors(ARGUMENTS.verifier, self.binary_path, self.text_path)

    def test_rangular_outputs_and_stage_timing(self) -> None:
        self.convert()
        for interaction_name, interactive_input in [
            ("full", "y\n"),
            ("first_order", "n\n1\n1\n1\n"),
        ]:
            generated_directories = []
            for input_mode in ("text", "descriptors"):
                for timing_enabled in (False, True):
                    run_directory = (
                        self.root / f"{interaction_name}_{input_mode}_{timing_enabled}"
                    )
                    run_directory.mkdir()
                    if input_mode == "text":
                        shutil.copyfile(self.text_path, run_directory / "rcsf.inp")
                    command = [str(ARGUMENTS.rangular)]
                    if input_mode == "descriptors":
                        command += ["--csf-descriptors", str(self.binary_path)]
                    if timing_enabled:
                        command += ["--timings", "timings.csv"]
                    result = subprocess.run(
                        command,
                        cwd=run_directory,
                        input=interactive_input,
                        text=True,
                        capture_output=True,
                        timeout=60,
                        check=False,
                        env={
                            **os.environ,
                            "OMP_NUM_THREADS": "1",
                            "OPENBLAS_NUM_THREADS": "1",
                        },
                    )
                    self.assertEqual(
                        result.returncode, 0, result.stdout + result.stderr
                    )
                    # Supply provenance for the independent comparator only after the binary-only run.
                    if input_mode == "descriptors":
                        self.assertFalse((run_directory / "rcsf.inp").exists())
                        shutil.copyfile(self.text_path, run_directory / "rcsf.inp")
                    outputs = {
                        path.name: path.read_bytes()
                        for path in run_directory.glob("mcp.*")
                    }
                    self.assertEqual(len(outputs), 6)
                    if generated_directories:
                        self.assertEqual(outputs, generated_directories[0][1])
                    generated_directories.append((run_directory, outputs))
                    if timing_enabled:
                        with (run_directory / "timings.csv").open() as stream:
                            timings = list(csv.DictReader(stream))
                        self.assertTrue(
                            all(float(item["wall_seconds"]) >= 0 for item in timings)
                        )
                        self.assertEqual(
                            sum(item["stage"] == "csf_load" for item in timings), 3
                        )
                        self.assertEqual(
                            sum(item["stage"] == "sort_and_output" for item in timings),
                            15,
                        )
                        self.assertEqual(
                            timings[-1]["stage"],
                            "program_total_including_interaction_input",
                        )
            # Independent MCP parser checks structure and nonzero coefficient coverage.
            compare_script = (
                Path(__file__).resolve().parents[1] / "scripts" / "mcp_compare.py"
            )
            comparison = subprocess.run(
                [
                    sys.executable,
                    str(compare_script),
                    "compare",
                    str(generated_directories[0][0]),
                    str(generated_directories[-1][0]),
                ],
                text=True,
                capture_output=True,
                timeout=30,
                check=False,
            )
            self.assertEqual(
                comparison.returncode, 0, comparison.stdout + comparison.stderr
            )
            report = json.loads(comparison.stdout)
            self.assertEqual(report["status"], "equal")
            self.assertTrue(
                any(channel["compared"] > 0 for channel in report["channels"]), report
            )


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
