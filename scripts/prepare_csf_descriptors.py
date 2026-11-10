#!/usr/bin/env python3
"""Convert bound rCSFs physical-seniority V2 Parquet to the suite stream format.

This is a container conversion. Fortran owns physical-state/coupling validation.
No rCSFs import, Rust build, or CSF text restoration is involved.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import struct
import sys
import time
import tomllib
from pathlib import Path

import numpy
import pyarrow as arrow
from pyarrow import parquet

MAGIC = b"GCSDESC1"
FORMAT_VERSION = 1
ENDIAN_MARKER = 0x01020304
PREFIX = struct.Struct("<8s6IQI64s")
ORBITAL = struct.Struct("<ii")
BLOCK = struct.Struct("<Qii")
MAX_COUNT = 2**31 - 1  # Current suite consumers use default INTEGER counts.
ORBITAL_LETTERS = "spdfghiklm"


def parse_orbital(label: str) -> tuple[int, int]:
    match = re.fullmatch(r"([1-9][0-9]*)([spdfghiklm])(-?)", label)
    if match is None:
        raise ValueError(f"invalid ordinary-GRASP orbital label: {label!r}")
    principal_number = int(match[1])
    orbital_momentum = ORBITAL_LETTERS.index(match[2])
    if principal_number <= orbital_momentum or principal_number > MAX_COUNT:
        raise ValueError(f"invalid principal quantum number: {label!r}")
    if match[3] and orbital_momentum == 0:
        raise ValueError("s- is not a relativistic subshell")
    relativistic_kappa = orbital_momentum if match[3] else -orbital_momentum - 1
    return principal_number, relativistic_kappa


def column_names(peel_count: int) -> list[str]:
    return [
        f"sub{shell_index}_{channel}"
        for shell_index in range(peel_count)
        for channel in ("n", "2j", "v", "2k")
    ] + ["total_two_j", "parity"]


def prepare_descriptors(
    descriptor_path: Path,
    header_path: Path,
    output_path: Path,
    *,
    batch_rows: int = 8192,
) -> dict[str, object]:
    start_time = time.perf_counter()
    if batch_rows <= 0:
        raise ValueError("batch_rows must be positive")
    header_bytes = header_path.read_bytes()
    header_digest = hashlib.sha256(header_bytes).hexdigest()
    header = tomllib.loads(header_bytes.decode("utf-8"))
    header_lines = header["header_info"]["header_lines"]
    if (
        not isinstance(header_lines, list)
        or len(header_lines) != 5
        or not all(isinstance(line, str) for line in header_lines)
    ):
        raise ValueError("header_info.header_lines must contain five strings")
    if [header_lines[index].strip() for index in (0, 2, 4)] != [
        "Core subshells:",
        "Peel subshells:",
        "CSF(s):",
    ]:
        raise ValueError("unexpected ordinary-CSF header identifiers")
    core_labels = header_lines[1].replace(",", " ").split()
    peel_labels = header_lines[3].replace(",", " ").split()
    orbitals = [parse_orbital(label) for label in core_labels + peel_labels]
    if not orbitals or len(set(orbitals)) != len(orbitals):
        raise ValueError("orbital table must be nonempty and contain no duplicates")
    if len(orbitals) > MAX_COUNT:
        raise ValueError("orbital count exceeds supported integer range")
    block_lengths = header["block_info"]["block_lengths"]
    if (
        not isinstance(block_lengths, list)
        or not block_lengths
        or any(
            type(count) is not int or not 0 < count <= MAX_COUNT
            for count in block_lengths
        )
    ):
        raise ValueError(
            "block_lengths must be a nonempty list of positive supported integers"
        )
    total_rows = sum(block_lengths)
    if total_rows > MAX_COUNT or len(block_lengths) > MAX_COUNT:
        raise ValueError("CSF/block count exceeds current suite integer range")
    peel_count = len(peel_labels)
    names = column_names(peel_count)
    source = parquet.ParquetFile(descriptor_path)
    # Rust writes the contract in the Parquet footer; ARROW:schema need not
    # carry it. Reading only schema_arrow.metadata misses genuine rCSFs files.
    metadata = source.metadata.metadata or {}
    required = {
        "descriptor_version": "2",
        "descriptor_dtype": "int8",
        "max_two_j": "126",
        "seniority_encoding": "physical",
        "channels_per_subshell": "4",
        "subshell_count": str(peel_count),
        "missing_sentinel": "-1",
        "normalized": "false",
        "source_header_sha256": header_digest,
    }
    for key, expected in required.items():
        if metadata.get(key.encode()) != expected.encode():
            raise ValueError(f"descriptor metadata requires {key}={expected}")
    for key, expected in {
        "peel_subshells": peel_labels,
        "feature_columns": names[:-2],
        "global_columns": names[-2:],
    }.items():
        if json.loads(metadata.get(key.encode(), b"null")) != expected:
            raise ValueError(
                f"descriptor metadata {key} disagrees with the header/schema"
            )
    if source.schema_arrow.names != names or any(
        field.type != arrow.int8() or field.nullable for field in source.schema_arrow
    ):
        raise ValueError("expected ordered, non-nullable Int8 V2 columns")
    if source.metadata.num_rows != total_rows:
        raise ValueError("block_lengths do not cover exactly the descriptor rows")

    shell_capacities = numpy.array(
        [2 * abs(kappa) for _, kappa in orbitals[len(core_labels) :]], dtype=numpy.int64
    )
    core_electrons = sum(2 * abs(kappa) for _, kappa in orbitals[: len(core_labels)])
    block_symmetries: list[tuple[int, int] | None] = [None] * len(block_lengths)
    block_index = 0
    rows_in_block = 0
    electron_count: int | None = None
    rows_written = 0
    # Exclusive creation prevents replacing either input or any previous result.
    with output_path.open("xb") as output:
        try:
            output.write(
                PREFIX.pack(
                    MAGIC,
                    FORMAT_VERSION,
                    ENDIAN_MARKER,
                    len(core_labels),
                    peel_count,
                    len(block_lengths),
                    len(names),
                    total_rows,
                    0,
                    header_digest.encode("ascii"),
                )
            )
            for principal_number, relativistic_kappa in orbitals:
                output.write(ORBITAL.pack(principal_number, relativistic_kappa))
            block_table_offset = output.tell()
            for row_count in block_lengths:
                output.write(BLOCK.pack(row_count, 0, 0))
            for batch in source.iter_batches(batch_size=batch_rows, columns=names):
                if batch.num_rows == 0:
                    continue
                if any(column.null_count for column in batch.columns):
                    raise ValueError("null descriptor value")
                values = numpy.column_stack(
                    [column.to_numpy() for column in batch.columns]
                )
                occupations = values[:, :-2:4]
                seniorities = values[:, 2:-2:4]
                angular_values = values[:, 1:-2:2]
                if (
                    numpy.any(occupations < 0)
                    or numpy.any(occupations > shell_capacities)
                    or numpy.any(seniorities < 0)
                    or numpy.any(angular_values < -1)
                    or numpy.any(angular_values > 126)
                    or numpy.any(values[:, -2] < 0)
                    or numpy.any(values[:, -2] > 126)
                    or numpy.any(numpy.abs(values[:, -1].astype(numpy.int16)) != 1)
                ):
                    raise ValueError(
                        f"invalid descriptor channel near row {rows_written + 1}"
                    )
                row_electrons = (
                    occupations.sum(axis=1, dtype=numpy.int64) + core_electrons
                )
                if electron_count is None:
                    electron_count = int(row_electrons[0])
                if not 0 < electron_count <= MAX_COUNT or numpy.any(
                    row_electrons != electron_count
                ):
                    raise ValueError("inconsistent or unsupported electron count")
                batch_offset = 0
                while batch_offset < batch.num_rows:
                    portion_rows = min(
                        batch.num_rows - batch_offset,
                        block_lengths[block_index] - rows_in_block,
                    )
                    portion = values[batch_offset : batch_offset + portion_rows]
                    symmetry = (int(portion[0, -2]), int(portion[0, -1]))
                    previous = block_symmetries[block_index]
                    if (previous is not None and previous != symmetry) or numpy.any(
                        portion[:, -2:] != symmetry
                    ):
                        raise ValueError(
                            f"mixed J/parity inside block {block_index + 1}"
                        )
                    block_symmetries[block_index] = symmetry
                    rows_in_block += portion_rows
                    batch_offset += portion_rows
                    if rows_in_block == block_lengths[block_index]:
                        block_index += 1
                        rows_in_block = 0
                output.write(values.tobytes(order="C"))
                rows_written += batch.num_rows
            if rows_written != total_rows or electron_count is None:
                raise ValueError("incomplete descriptor data")
            output.seek(40)  # electron count in the version-1 prefix
            output.write(struct.pack("<I", electron_count))
            output.seek(block_table_offset)
            for row_count, symmetry in zip(
                block_lengths, block_symmetries, strict=True
            ):
                if symmetry is None:
                    raise ValueError("missing block symmetry")
                output.write(BLOCK.pack(row_count, *symmetry))
        except BaseException:
            output.close()
            output_path.unlink()
            raise
    return {
        "format_version": FORMAT_VERSION,
        "rows": total_rows,
        "blocks": len(block_lengths),
        "orbitals": len(orbitals),
        "header_sha256": header_digest,
        "output_bytes": output_path.stat().st_size,
        "preprocess_wall_seconds": time.perf_counter() - start_time,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("descriptors", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--header", type=Path, required=True)
    parser.add_argument("--batch-rows", type=int, default=8192)
    arguments = parser.parse_args()
    try:
        report = prepare_descriptors(
            arguments.descriptors,
            arguments.header,
            arguments.output,
            batch_rows=arguments.batch_rows,
        )
    except (OSError, ValueError, KeyError, TypeError, arrow.ArrowException) as error:
        print(f"CSF descriptor preparation failed: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
