#!/usr/bin/env python3
"""Validate, inspect, export and compare sorted GRASP rangular MCP files.

NumPy handles coefficient arrays in batches. No GRASP executable is needed.
See docs/implemented/mcp_compare.md for the format and comparison contract.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import re
import struct
import sys
import time
from collections.abc import Iterator
from contextlib import ExitStack
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
from numpy.typing import NDArray

Array = NDArray[Any]
CHUNK = 1_000_000


class MCPError(ValueError):
    """Invalid, incomplete or unsupported MCP data."""


class Records:
    """Sequential unformatted records, including GNU continuation subrecords."""

    def __init__(self, path: Path, stream: BinaryIO, structure: bool):
        self.path, self.stream = path, stream
        self.size = path.stat().st_size
        prefix = stream.read(64)
        stream.seek(0)
        candidates = []
        for endian in ("<", ">"):
            for width, code in ((4, "i"), (8, "q")):
                if len(prefix) < width:
                    continue
                n = struct.unpack_from(endian + code, prefix)[0]
                allowed = (12, 24) if structure else (15, 27)
                if n not in allowed or len(prefix) < 2 * width + n:
                    continue
                if struct.unpack_from(endian + code, prefix, width + n)[0] == n and (
                    structure or prefix[width : width + 3] == b"MCP"
                ):
                    candidates.append((endian, width, code, n))
        if len(candidates) != 1:
            self.fail("cannot detect record markers / endian / integer width")
        self.endian, self.width, code, n = candidates[0]
        self.marker = struct.Struct(self.endian + code)
        self.int_width = n // 3 if structure else (n - 3) // 3
        self.idtype = np.dtype(f"{self.endian}i{self.int_width}")

    def fail(self, message: str) -> None:
        raise MCPError(f"{self.path} at byte {self.stream.tell()}: {message}")

    def record(self, expected: int | None = None) -> bytes:
        parts = []
        total = 0
        previous = False
        while True:
            raw = self.stream.read(self.width)
            if len(raw) != self.width:
                self.fail("truncated/missing record marker")
            leading = self.marker.unpack(raw)[0]
            size = abs(leading)
            if size > self.size - self.stream.tell() - self.width:
                self.fail("truncated record payload")
            if expected is not None and total + size > expected:
                self.fail(f"record exceeds expected {expected} bytes")
            payload = self.stream.read(size)
            raw = self.stream.read(self.width)
            if len(payload) != size or len(raw) != self.width:
                self.fail("truncated record")
            trailing = self.marker.unpack(raw)[0]
            if abs(trailing) != size or (trailing < 0) != previous:
                self.fail("mismatched record markers / continuation signs")
            parts.append(payload)
            total += size
            if leading >= 0:
                break
            previous = True
        if expected is not None and total != expected:
            self.fail(f"record has {total} bytes, expected {expected}")
        return parts[0] if len(parts) == 1 else b"".join(parts)

    def ints(self, count: int) -> tuple[int, ...]:
        return tuple(
            map(int, np.frombuffer(self.record(count * self.int_width), self.idtype))
        )

    def tagged(self, count: int) -> tuple[int, ...]:
        data = self.record(3 + count * self.int_width)
        if data[:3] != b"MCP":
            self.fail("missing MCP signature (temporary/CSFG formats are unsupported)")
        return tuple(map(int, np.frombuffer(data, self.idtype, offset=3)))

    def eof(self) -> None:
        if self.stream.read(1):
            self.fail("unexpected trailing data")


@dataclass(frozen=True)
class Header:
    nblock: int
    rank: int
    nprocs: int
    nelec: int
    ncf: int
    nw: int
    diag: bool
    iccut: int
    lfordr: bool

    def semantics(self) -> tuple[Any, ...]:
        return (
            self.nblock,
            self.nelec,
            self.ncf,
            self.nw,
            self.diag,
            self.iccut,
            self.lfordr,
        )


def header(reader: Records) -> Header:
    nb, rank, procs = reader.tagged(3)
    nelec, ncf, nw = reader.ints(3)
    diag, iccut, lfordr = reader.ints(3)
    if not (
        nb > 0 and 0 <= rank < procs and nelec > 0 and ncf > 0 and nw > 0 and iccut >= 0
    ):
        reader.fail("invalid global header values")
    if diag not in (0, 1, -1) or lfordr not in (0, 1, -1):
        reader.fail("unsupported LOGICAL representation")
    return Header(nb, rank, procs, nelec, ncf, nw, bool(diag), iccut, bool(lfordr))


@dataclass
class SparseBlock:
    ncf: int
    rank: int
    nprocs: int
    ends: Array
    rows: Array

    def columns(self, positions: Array) -> Array:
        return (
            np.searchsorted(self.ends, positions, side="right") * self.nprocs
            + self.rank
            + 1
        )

    def pairs(self) -> Array:
        result = np.empty(self.rows.size, dtype=np.int64)
        for start in range(0, self.rows.size, CHUNK):
            stop = min(start + CHUNK, self.rows.size)
            cols = self.columns(np.arange(start, stop, dtype=np.int64))
            result[start:stop] = cols * (self.ncf + 1) + self.rows[start:stop]
        return result


@dataclass
class Shard:
    prefix: Path
    header: Header
    blocks: list[SparseBlock]


@dataclass(frozen=True)
class IndexedChannel:
    first: int
    groups: int
    contributions: int
    position: int


@dataclass
class IndexedShard(Shard):
    directory_position: int
    channels: list[list[IndexedChannel]]
    endian: str
    radix: int
    generation: int
    basis_id: bytes
    cutoffs: tuple[int, ...]


def read_indexed_shard(path: Path) -> tuple[IndexedShard, tuple[Any, ...]]:
    """Validate the stream envelope and load the compact block directory."""
    with path.open("rb") as stream:
        size = path.stat().st_size
        def fail(message: str) -> None:
            raise MCPError(f"{path} at byte {stream.tell()}: {message}")

        def unpack(fields: str) -> tuple[Any, ...]:
            layout = struct.Struct(endian + fields)
            data = stream.read(layout.size)
            if len(data) != layout.size:
                fail("truncated indexed MCP")
            return layout.unpack(data)

        prefix = stream.read(16)
        if prefix[:8] != b"GMCPIDX1" or prefix[12:16] not in (b"\x04\x03\x02\x01", b"\x01\x02\x03\x04"):
            fail("invalid indexed MCP signature or byte order")
        endian = "<" if prefix[12:16] == b"\x04\x03\x02\x01" else ">"
        stream.seek(8)
        version, _, i32, i64, real, radix, integer_format, real_format = unpack("8i")
        if (version, i32, i64, real, integer_format, real_format) != (1, 4, 8, 8, 1, 1) or radix < 2:
            fail("unsupported indexed MCP representation")
        rank, ranks, core, blocks, kmax, electrons, csfs, nw, flags, generation, directory, groups = unpack("12q")
        basis_id, = unpack("64s")
        if not (0 <= rank < ranks and 0 < nw < radix and 0 <= core <= nw and 1 <= blocks <= 50
                and 0 <= kmax <= 67 and electrons > 0 and 0 <= flags <= 3
                and groups >= 0 and 200 < directory <= size - 31
                and directory - 1 + 32 * (groups + 1) == size):
            fail("invalid header dimensions or truncated indexed MCP directory")
        if not re.fullmatch(rb"[a-f0-9]{64}", basis_id):
            fail("invalid basis identifier")
        expected_position = 201 + blocks * (48 + 24 * (kmax + 2))
        expected_group = 1
        channels, sparse_data, ids, ncfs, cutoffs = [], [], [], [], []
        for _ in range(blocks):
            tag, ncf, cutoff, position, ncols, nrows = unpack("8s5q")
            columns_owned = max(0, (ncf - rank - 1) // ranks + 1)
            if not (0 < ncf < 2**31 and 0 <= cutoff <= ncf and position == expected_position
                    and ncols == columns_owned and ncols <= nrows < 2**31
                    and position + 4 * (nrows + ncols) <= directory):
                fail("invalid indexed block structure")
            ncfs.append(ncf)
            ids.append(tag.decode("ascii", errors="replace").strip())
            cutoffs.append(cutoff)
            sparse_data.append((position, ncols, nrows))
            expected_position += 4 * (ncols + nrows)
            first, counts, totals = [unpack(f"{kmax + 2}q") for _ in range(3)]
            block_channels = []
            for first_group, count, total in zip(first, counts, totals, strict=True):
                if not (first_group == expected_group and 0 <= count <= groups - expected_group + 1
                        and count <= total and expected_position + 16 * total <= directory):
                    fail("invalid indexed channel extent")
                block_channels.append(IndexedChannel(first_group, count, total, expected_position))
                expected_group += count
                expected_position += 16 * total
            channels.append(block_channels)
        if expected_group != groups + 1 or expected_position != directory or sum(ncfs) != csfs:
            fail("indexed directory totals do not match header")
        stream.seek(size - 32)
        if unpack("8s3q") != (b"MCPDONE1", generation, directory, groups):
            fail("incomplete indexed MCP shard")
        dtype = np.dtype(endian + "i4")
        sparse = []
        for ncf, (position, ncols, nrows) in zip(ncfs, sparse_data, strict=True):
            stream.seek(position - 1)
            data = np.fromfile(stream, dtype, ncols + nrows)
            if data.size != ncols + nrows:
                fail("truncated indexed sparse structure")
            ends, rows = data[:ncols], data[ncols:]
            if ncols:
                cols = np.arange(rank + 1, ncf + 1, ranks)
                if ends[0] < 1 or ends[-1] != nrows or np.any(ends[1:] <= ends[:-1]):
                    fail("invalid indexed sparse column endpoints")
                if np.any(rows[ends - 1] != cols):
                    fail("missing diagonal at end of indexed column")
            elif nrows:
                fail("indexed sparse rows on an empty rank")
            block = SparseBlock(ncf, rank, ranks, ends, rows)
            pairs = block.pairs()
            for start in range(0, nrows, CHUNK):
                stop = min(start + CHUNK, nrows)
                cols = block.columns(np.arange(start, stop))
                if np.any(rows[start:stop] < 1) or np.any(rows[start:stop] > cols):
                    fail("indexed row outside lower CSF triangle")
            pairs.sort()
            if np.any(pairs[1:] == pairs[:-1]):
                fail("duplicate indexed sparse CSF pair")
            sparse.append(block)
    h = Header(blocks, rank, ranks, electrons, csfs, nw, bool(flags & 1), cutoffs[0], bool(flags & 2))
    return IndexedShard(path, h, sparse, directory, channels, endian, radix, generation, basis_id, tuple(cutoffs)), (
        core, tuple(ncfs), tuple(ids), kmax, h.semantics()
    )


def indexed_coefficient_groups(
    stream: BinaryIO, shard: IndexedShard, unit: int, block: int, radix: int
) -> Iterator[tuple[int, Array, Array]]:
    if radix != shard.radix:
        raise MCPError(f"{shard.prefix}: orbital label radix is {shard.radix}; check --keyorb")
    channel = shard.channels[block-1][unit-31]
    sparse = shard.blocks[block-1]
    dtype = np.dtype([("column", shard.endian+"i4"), ("index", shard.endian+"i4"), ("value", shard.endian+"f8")])
    remaining, position, previous = channel.contributions, channel.position, 0
    for index in range(channel.first, channel.first + channel.groups):
        stream.seek(shard.directory_position - 1 + (index - 1)*32)
        payload = stream.read(32)
        if len(payload) != 32:
            raise MCPError(f"{shard.prefix}: truncated indexed group directory")
        offset, length, label, count = struct.unpack(shard.endian+"4q", payload)
        if offset != position or not 0 < count <= remaining or length != 16*count or label <= previous:
            raise MCPError(f"{shard.prefix}: invalid indexed integral group")
        if any(i < 1 or i > shard.header.nw for i in orbitals(label, unit, radix)):
            raise MCPError(f"{shard.prefix}: indexed group label outside orbital basis")
        stream.seek(offset - 1)
        data = np.fromfile(stream, dtype, count)
        if data.size != count:
            raise MCPError(f"{shard.prefix}: truncated indexed coefficient payload")
        pairs, values = checked_contributions(sparse, data, shard.prefix)
        yield label, pairs, values
        remaining -= count
        position += length
        previous = label
    if remaining:
        raise MCPError(f"{shard.prefix}: indexed channel contribution total mismatch")


class Dataset:
    """A complete serial set or a complete collection of MPI rank files."""

    def __init__(self, directory: Path):
        self.directory = directory
        channels = {
            p.relative_to(directory)
            for pattern in ("mcp*.*", "[0-9][0-9][0-9]/mcp.*")
            for p in directory.glob(pattern)
            if re.fullmatch(r"mcp(?:\d+)?\.\d+", p.name)
        }
        indexed_paths = sorted(
            path for pattern in ("mcp*.mcpidx", "[0-9][0-9][0-9]/mcp.mcpidx") for path in directory.glob(pattern)
            if re.fullmatch(r"mcp(?:[0-9]+)?\.mcpidx", path.name)
        )
        if indexed_paths:
            if channels:
                raise MCPError(f"{directory}: mixed indexed and sequential MCP file sets")
            self.load_indexed(indexed_paths)
            return
        files = [directory / p for p in sorted(channels) if p.suffix == ".30"]
        if not files:
            raise MCPError(f"{directory}: no mcp.30 or mcpNNN.30 files")
        if (directory / "mcp.30") in files and len(files) != 1:
            raise MCPError(
                f"{directory}: mixed serial and MPI files; use separate directories"
            )
        self.shards: list[Shard] = []
        self.paths: list[Path] = []
        for path in files:
            with path.open("rb") as stream:
                r = Records(path, stream, True)
                ncore, nblock, kmax = r.ints(3)
                if not (
                    ncore >= 0
                    and 1 <= nblock <= r.size // r.int_width
                    and 0 <= kmax <= 67
                ):
                    r.fail("invalid structure preamble")
                ncfs = r.ints(nblock)
                try:
                    raw_ids = r.record(8 * nblock).decode("ascii")
                except UnicodeDecodeError as exc:
                    raise MCPError(f"{path}: non-ASCII block IDs") from exc
                ids = tuple(
                    raw_ids[i : i + 8].strip() for i in range(0, len(raw_ids), 8)
                )
                h = header(r)
                if (
                    h.nblock != nblock
                    or sum(ncfs) != h.ncf
                    or min(ncfs) <= 0
                    or ncore > h.nw
                ):
                    r.fail("inconsistent block counts / core size")
                if max(ncfs) > 3_037_000_498:
                    r.fail("CSF count exceeds safe int64 pair-key capacity")
                meta = (ncore, ncfs, ids, kmax, h.semantics())
                if not self.shards:
                    self.metadata = meta
                    self.ncore, self.ncfs, self.ids, self.kmax = ncore, ncfs, ids, kmax
                    self.header = h
                elif meta != self.metadata or h.nprocs != self.header.nprocs:
                    r.fail("inconsistent headers between MPI ranks")
                if path.stem != "mcp" and int(path.stem[3:]) != h.rank:
                    r.fail("filename rank does not match header")
                if path.parent != directory and int(path.parent.name) != h.rank:
                    r.fail("rank directory does not match header")
                blocks = []
                for block, ncf in enumerate(ncfs, 1):
                    if r.tagged(2) != (block, ncf):
                        r.fail(f"unexpected structure block {block}")
                    (nnz,) = r.ints(1)
                    ncolumns = max(0, (ncf - h.rank - 1) // h.nprocs + 1)
                    if nnz < ncolumns:
                        r.fail("too few sparse entries (diagonals must be present)")
                    data = np.frombuffer(
                        r.record((ncolumns + nnz) * r.int_width), r.idtype
                    )
                    columns = np.arange(h.rank + 1, ncf + 1, h.nprocs, dtype=np.int64)
                    ends, rows = data[: columns.size], data[columns.size :]
                    if columns.size:
                        if (
                            ends[0] <= 0
                            or ends[-1] != nnz
                            or np.any(ends[1:] <= ends[:-1])
                        ):
                            r.fail("invalid sparse column endpoints")
                        if np.any(rows[ends - 1] != columns):
                            r.fail("missing diagonal at end of sparse column")
                    elif nnz:
                        r.fail("nonempty structure on rank with no columns")
                    sparse = SparseBlock(ncf, h.rank, h.nprocs, ends, rows)
                    for start in range(0, nnz, CHUNK):
                        stop = min(start + CHUNK, nnz)
                        cols = sparse.columns(np.arange(start, stop, dtype=np.int64))
                        if np.any(rows[start:stop] < 1) or np.any(
                            rows[start:stop] > cols
                        ):
                            r.fail("row outside lower CSF triangle")
                    pairs = sparse.pairs()
                    pairs.sort()
                    if np.any(pairs[1:] == pairs[:-1]):
                        r.fail("duplicate sparse CSF pair")
                    blocks.append(sparse)
                r.eof()
                self.shards.append(Shard(path.with_suffix(""), h, blocks))
        self.shards.sort(key=lambda s: s.header.rank)
        if [s.header.rank for s in self.shards] != list(range(self.header.nprocs)):
            raise MCPError(f"{directory}: missing or duplicate MPI ranks")
        expected = {
            s.prefix.with_suffix(f".{unit}").relative_to(directory)
            for s in self.shards
            for unit in range(30, self.kmax + 33)
        }
        if channels != expected:
            raise MCPError(
                f"{directory}: incomplete/unexpected file set; missing={sorted(expected - channels)}, extra={sorted(channels - expected)}"
            )
        self.paths = [directory / name for name in sorted(expected)]

    def load_indexed(self, paths: list[Path]) -> None:
        self.shards = []
        first, self.metadata = read_indexed_shard(paths[0])
        self.ncore, self.ncfs, self.ids, self.kmax, _ = self.metadata
        self.header = first.header
        self.paths = paths
        for path in paths:
            shard, metadata = (first, self.metadata) if path == paths[0] else read_indexed_shard(path)
            if (metadata != self.metadata or shard.header.nprocs != self.header.nprocs
                    or (shard.generation, shard.basis_id, shard.cutoffs, shard.radix, shard.endian)
                    != (first.generation, first.basis_id, first.cutoffs, first.radix, first.endian)):
                raise MCPError(f"{path}: inconsistent indexed shard metadata")
            if path.parent != self.directory and int(path.parent.name) != shard.header.rank:
                raise MCPError(f"{path}: rank directory does not match header")
            if path.stem != "mcp" and int(path.stem[3:]) != shard.header.rank:
                raise MCPError(f"{path}: filename rank does not match header")
            self.shards.append(shard)
        self.shards.sort(key=lambda shard: shard.header.rank)
        if len(self.shards) != self.header.nprocs or any(
            shard.header.rank != index for index, shard in enumerate(self.shards)
        ):
            raise MCPError(f"{self.directory}: missing or duplicate indexed ranks")

    def summary(self) -> dict[str, Any]:
        return {
            "directory": str(self.directory),
            "ncore": self.ncore,
            "header": asdict(self.header),
            "kmax": self.kmax,
            "blocks": [
                {
                    "block": b + 1,
                    "id": ident,
                    "ncf": ncf,
                    "sparse_entries": sum(s.blocks[b].rows.size for s in self.shards),
                }
                for b, (ident, ncf) in enumerate(zip(self.ids, self.ncfs))
            ],
            "files": len(self.paths),
            "bytes": sum(p.stat().st_size for p in self.paths),
        }

    def structure(self, block: int) -> Array:
        pairs = np.concatenate([s.blocks[block - 1].pairs() for s in self.shards])
        pairs.sort()
        return pairs

    def groups(self, unit: int, keyorb: int) -> Iterator[tuple[int, int, Array, Array]]:
        """Yield (block, label, CSF pair keys, coefficients), one label at a time."""
        with ExitStack() as stack:
            readers: list[Records] = []
            indexed = isinstance(self.shards[0], IndexedShard)
            files: list[BinaryIO] = []
            for shard in self.shards:
                if indexed:
                    files.append(stack.enter_context(shard.prefix.open("rb")))
                else:
                    path = shard.prefix.with_suffix(f".{unit}")
                    r = Records(path, stack.enter_context(path.open("rb")), False)
                    if header(r) != shard.header:
                        r.fail("coefficient header differs from mcp.30")
                    readers.append(r)
            for block in range(1, self.header.nblock + 1):
                if indexed:
                    streams = [
                        iter(indexed_coefficient_groups(file, shard, unit, block, keyorb))
                        for file, shard in zip(files, self.shards, strict=True)
                        if isinstance(shard, IndexedShard)
                    ]
                else:
                    streams = [
                        iter(coefficient_groups(r, s.blocks[block - 1], unit, block, keyorb, self.header.nw))
                        for r, s in zip(readers, self.shards, strict=True)
                    ]
                heap = []
                for i, source in enumerate(streams):
                    item = next(source, None)
                    if item is not None:
                        heapq.heappush(heap, (item[0], i, item[1], item[2]))
                while heap:
                    label = heap[0][0]
                    pairs, values = [], []
                    while heap and heap[0][0] == label:
                        _, i, p, v = heapq.heappop(heap)
                        pairs.append(p)
                        values.append(v)
                        item = next(streams[i], None)
                        if item is not None:
                            heapq.heappush(heap, (item[0], i, item[1], item[2]))
                    yield block, label, np.concatenate(pairs), np.concatenate(values)
            for r in readers:
                r.eof()


def orbitals(label: int, unit: int, keyorb: int) -> tuple[int, ...]:
    if unit == 31:
        b, a = divmod(label, keyorb)
        return a, b
    ac, bd = divmod(label, keyorb * keyorb)
    a, c = divmod(ac, keyorb)
    b, d = divmod(bd, keyorb)
    return a, b, c, d


def coefficient_groups(
    r: Records, sparse: SparseBlock, unit: int, block: int, keyorb: int, nw: int
) -> Iterator[tuple[int, Array, Array]]:
    b, ncf, count = r.tagged(3)
    if (b, ncf) != (block, sparse.ncf) or count < 0:
        r.fail(f"invalid coefficient block {block}")
    dtype = np.dtype(
        [("column", r.idtype), ("index", r.idtype), ("value", r.endian + "f8")]
    )
    total, previous = 0, 0
    while True:
        label, size = r.ints(2)
        if (label, size) == (0, 0):
            break
        if label <= previous or size <= 0 or size > count - total:
            r.fail("invalid label order / contribution count")
        if any(i < 1 or i > nw for i in orbitals(label, unit, keyorb)):
            r.fail(f"label {label} has invalid orbital indices; check --keyorb")
        previous = label
        data = np.frombuffer(r.record(size * dtype.itemsize), dtype)
        pairs, values = checked_contributions(sparse, data, r.path)
        total += size
        yield label, pairs, values
    if total != count:
        r.fail(f"NCOEFF={count}, parsed {total} contributions")


def checked_contributions(sparse: SparseBlock, data: Array, path: Path) -> tuple[Array, Array]:
    positions = data["index"].astype(np.int64) - 1
    if np.any(positions < 0) or np.any(positions >= sparse.rows.size):
        raise MCPError(f"{path}: coefficient index outside sparse structure")
    cols = sparse.columns(positions)
    if np.any(cols != data["column"]):
        raise MCPError(f"{path}: coefficient column disagrees with sparse index")
    if not np.all(np.isfinite(data["value"])):
        raise MCPError(f"{path}: non-finite angular coefficient")
    return cols * (sparse.ncf + 1) + sparse.rows[positions], data["value"]


def canonical(pairs: Array, values: Array, sum_duplicates: bool) -> tuple[Array, Array]:
    order = np.lexsort((values, pairs))
    pairs, values = pairs[order], values[order]
    first = (
        np.r_[0, np.flatnonzero(pairs[1:] != pairs[:-1]) + 1]
        if pairs.size
        else np.array([], dtype=np.int64)
    )
    keys = np.empty(
        first.size if sum_duplicates else pairs.size,
        dtype=[("pair", "i8"), ("occurrence", "i8")],
    )
    if sum_duplicates:
        keys["pair"] = pairs[first]
        keys["occurrence"] = 0
        with np.errstate(over="ignore", invalid="ignore"):
            values = np.add.reduceat(values, first) if first.size else values
        if not np.all(np.isfinite(values)):
            raise MCPError("non-finite sum of duplicate coefficients")
    else:
        keys["pair"] = pairs
        keys["occurrence"] = np.arange(pairs.size) - np.repeat(
            first, np.diff(np.r_[first, pairs.size])
        )
    return keys, values


def location(
    unit: int, block: int, label: int, pair: int, ncf: int, keyorb: int
) -> dict[str, Any]:
    col, row = divmod(pair, ncf + 1)
    return {
        "unit": unit,
        "kind": "T" if unit == 31 else "V",
        "k": None if unit == 31 else unit - 32,
        "block": block,
        "label": label,
        "orbitals": orbitals(label, unit, keyorb),
        "column_csf": col,
        "row_csf": row,
    }


def compare(
    reference: Dataset,
    candidate: Dataset,
    *,
    atol: float = 0.0,
    rtol: float = 0.0,
    keyorb: int = 215,
    sum_duplicates: bool = False,
    max_diffs: int = 20,
) -> dict[str, Any]:
    if not all(math.isfinite(v) and v >= 0 for v in (atol, rtol)) or max_diffs < 0:
        raise MCPError("tolerances and max_diffs must be finite and nonnegative")
    report: dict[str, Any] = {
        "status": "equal",
        "reference": reference.summary(),
        "candidate": candidate.summary(),
        "atol": atol,
        "rtol": rtol,
        "keyorb": keyorb,
        "sum_duplicates": sum_duplicates,
        "differences": [],
        "difference_count": 0,
        "channels": [],
    }

    def difference(item: dict[str, Any], count: int = 1) -> None:
        report["status"] = "different"
        report["difference_count"] += count
        if len(report["differences"]) < max_diffs:
            report["differences"].append(item)

    indexed_reference, indexed_candidate = reference.shards[0], candidate.shards[0]
    if isinstance(indexed_reference, IndexedShard) and isinstance(indexed_candidate, IndexedShard):
        if indexed_reference.cutoffs != indexed_candidate.cutoffs:
            difference({"type": "interaction_cutoffs", "reference": indexed_reference.cutoffs, "candidate": indexed_candidate.cutoffs})
        if indexed_reference.basis_id != indexed_candidate.basis_id:
            raise MCPError("indexed MCP basis identifiers differ")
    if reference.metadata != candidate.metadata:
        difference(
            {
                "type": "metadata",
                "reference": reference.metadata,
                "candidate": candidate.metadata,
            }
        )
        # Fully parse both sides even on a metadata mismatch: malformed files are errors.
        for dataset in (reference, candidate):
            inspect(dataset, keyorb)
        return report
    for block, ncf in enumerate(reference.ncfs, 1):
        left, right = reference.structure(block), candidate.structure(block)
        for side, a, b in (
            ("reference_only", left, right),
            ("candidate_only", right, left),
        ):
            missing = np.setdiff1d(a, b, assume_unique=True)
            if missing.size:
                difference(
                    {
                        "type": "sparse_structure",
                        "side": side,
                        "block": block,
                        "count": int(missing.size),
                        "first_pair": divmod(int(missing[0]), ncf + 1),
                    },
                    int(missing.size),
                )
    empty = np.array([], dtype=np.int64), np.array([], dtype=np.float64)
    for unit in range(31, reference.kmax + 33):
        stats: dict[str, Any] = {
            "unit": unit,
            "reference_contributions": 0,
            "candidate_contributions": 0,
            "compared": 0,
            "missing": 0,
            "mismatched": 0,
            "max_abs_error": 0.0,
            "max_rel_error": 0.0,
        }
        li, ri = (
            iter(reference.groups(unit, keyorb)),
            iter(candidate.groups(unit, keyorb)),
        )
        l, r = next(li, None), next(ri, None)
        while l is not None or r is not None:
            key = min(item[:2] for item in (l, r) if item is not None)
            lp, lv = l[2:] if l is not None and l[:2] == key else empty
            rp, rv = r[2:] if r is not None and r[:2] == key else empty
            stats["reference_contributions"] += len(lv)
            stats["candidate_contributions"] += len(rv)
            lk, lv = canonical(lp, lv, sum_duplicates)
            rk, rv = canonical(rp, rv, sum_duplicates)
            if np.array_equal(lk, rk):
                common, a, b = lk, lv, rv
            else:
                common, lix, rix = np.intersect1d(
                    lk, rk, assume_unique=True, return_indices=True
                )
                a, b = lv[lix], rv[rix]
                for side, keys, other in (
                    ("reference_only", lk, rk),
                    ("candidate_only", rk, lk),
                ):
                    missing = np.setdiff1d(keys, other, assume_unique=True)
                    if missing.size:
                        stats["missing"] += int(missing.size)
                        loc = location(
                            unit,
                            *key,
                            int(missing[0]["pair"]),
                            reference.ncfs[key[0] - 1],
                            keyorb,
                        )
                        difference(
                            {
                                "type": "contribution_missing",
                                "side": side,
                                "count": int(missing.size),
                                **loc,
                            },
                            int(missing.size),
                        )
            stats["compared"] += len(a)
            # longdouble has extra range on Linux. Use scaled arithmetic as a
            # fallback on platforms where it aliases float64 (notably Windows).
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                delta = np.abs(a.astype(np.longdouble) - b.astype(np.longdouble))
                scale = np.maximum(
                    np.abs(a.astype(np.longdouble)), np.abs(b.astype(np.longdouble))
                )
                relative = np.divide(
                    delta, scale, out=np.zeros_like(delta), where=scale != 0
                )
                overflow = ~np.isfinite(delta)
                relative[overflow] = (
                    np.abs(a[overflow]) / scale[overflow]
                    + np.abs(b[overflow]) / scale[overflow]
                )
                tolerance = np.longdouble(atol) + np.longdouble(rtol) * scale
                outside = delta > tolerance
                fallback = overflow | ~np.isfinite(tolerance)
                outside[fallback] = relative[fallback] > rtol + atol / scale[fallback]
            bad = np.flatnonzero(outside)
            if delta.size:
                stats["max_abs_error"] = max(stats["max_abs_error"], delta.max())
                stats["max_rel_error"] = max(
                    stats["max_rel_error"], float(relative.max())
                )
            if bad.size:
                stats["mismatched"] += int(bad.size)
                for i in bad[: max(1, max_diffs - len(report["differences"]))]:
                    loc = location(
                        unit,
                        *key,
                        int(common[i]["pair"]),
                        reference.ncfs[key[0] - 1],
                        keyorb,
                    )
                    difference(
                        {
                            "type": "coefficient",
                            **loc,
                            "occurrence": int(common[i]["occurrence"]) + 1,
                            "reference": float(a[i]),
                            "candidate": float(b[i]),
                            "abs_error": json_number(delta[i]),
                            "rel_error": float(relative[i]),
                        },
                        0,
                    )
                report["difference_count"] += int(bad.size)
            if l is not None and l[:2] == key:
                l = next(li, None)
            if r is not None and r[:2] == key:
                r = next(ri, None)
        stats["max_abs_error"] = json_number(stats["max_abs_error"])
        report["channels"].append(stats)
    return report


def json_number(value: Any) -> float | str:
    """Preserve out-of-float64-range errors without invalid JSON infinity."""
    return float(value) if abs(value) <= np.finfo(np.float64).max else str(value)


def inspect(dataset: Dataset, keyorb: int) -> dict[str, Any]:
    result = dataset.summary()
    result["status"] = "valid"
    result["channels"] = []
    for unit in range(31, dataset.kmax + 33):
        groups, count, max_abs = 0, 0, 0.0
        for _, _, _, values in dataset.groups(unit, keyorb):
            groups += 1
            count += len(values)
            max_abs = max(max_abs, float(np.max(np.abs(values))))
        result["channels"].append(
            {
                "unit": unit,
                "label_groups": groups,
                "contributions": count,
                "max_abs": max_abs,
            }
        )
    return result


def input_identity(reference: Path, candidate: Path, explicit: bool) -> dict[str, Any]:
    a, b = reference / "rcsf.inp", candidate / "rcsf.inp"
    if a.is_file() and b.is_file():
        with a.open("rb") as stream:
            ah = hashlib.file_digest(stream, "sha256").hexdigest()
        with b.open("rb") as stream:
            bh = hashlib.file_digest(stream, "sha256").hexdigest()
        if ah != bh:
            raise MCPError(
                "rcsf.inp differs: CSF/orbital ordering cannot be assumed identical"
            )
        return {"status": "sha256_identical", "sha256": ah}
    if not explicit:
        raise MCPError(
            "rcsf.inp missing on one/both sides; supply the inputs or use --assume-same-csfs after checking CSF/orbital order"
        )
    return {"status": "user_assumed_same_csfs"}


def nonnegative(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("must be finite and nonnegative")
    return number


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("inspect", "compare", "dump"):
        p = sub.add_parser(command)
        p.add_argument(
            "reference", type=Path, help="directory containing a complete MCP set"
        )
        p.add_argument(
            "--keyorb",
            type=int,
            default=215,
            help="label radix used by the producer (default: 215)",
        )
        if command == "compare":
            p.add_argument("candidate", type=Path)
            p.add_argument("--atol", type=nonnegative, default=0.0)
            p.add_argument("--rtol", type=nonnegative, default=0.0)
            p.add_argument("--sum-duplicates", action="store_true")
            p.add_argument("--assume-same-csfs", action="store_true")
            p.add_argument("--max-diffs", type=int, default=20)
        if command == "dump":
            p.add_argument("--unit", type=int, required=True)
            p.add_argument("--block", type=int)
            p.add_argument("--label", type=int)
            p.add_argument(
                "--limit", type=int, default=100, help="maximum JSONL rows (0: all)"
            )
        else:
            p.add_argument("--json", type=Path, help="also save the JSON report")
    args = parser.parse_args(argv)
    start = time.perf_counter()
    try:
        if args.keyorb < 2:
            raise MCPError("--keyorb must be >= 2")
        if args.command == "compare" and args.max_diffs < 0:
            raise MCPError("--max-diffs must be nonnegative")
        dataset = Dataset(args.reference)
        if dataset.header.nw >= args.keyorb:
            raise MCPError("NW must be less than --keyorb")
        other: Dataset | None = None
        if args.command == "compare":
            other = Dataset(args.candidate)
            if other.header.nw >= args.keyorb:
                raise MCPError("candidate NW must be less than --keyorb")
            first, second = dataset.shards[0], other.shards[0]
            if isinstance(first, IndexedShard) and isinstance(second, IndexedShard):
                if first.basis_id != second.basis_id:
                    raise MCPError("indexed MCP basis identifiers differ")
                identity = {"status": "indexed_basis_sha256_identical", "sha256": first.basis_id.decode("ascii")}
            else:
                identity = input_identity(args.reference,args.candidate,args.assume_same_csfs)
            result = compare(
                dataset,
                other,
                atol=args.atol,
                rtol=args.rtol,
                keyorb=args.keyorb,
                sum_duplicates=args.sum_duplicates,
                max_diffs=args.max_diffs,
            )
            result["csf_identity"] = identity
        elif args.command == "inspect":
            result = inspect(dataset, args.keyorb)
        else:
            if not 31 <= args.unit <= dataset.kmax + 32 or args.limit < 0:
                raise MCPError("invalid --unit or --limit")
            if args.block is not None and not 1 <= args.block <= dataset.header.nblock:
                raise MCPError("--block outside available blocks")
            emitted = 0
            # Consume the complete selected channel even after the output limit,
            # so truncation/invalid later blocks are not hidden by a preview.
            for block, label, pairs, values in dataset.groups(args.unit, args.keyorb):
                if args.block not in (None, block) or args.label not in (None, label):
                    continue
                remaining = (
                    len(values) if args.limit == 0 else max(0, args.limit - emitted)
                )
                for pair, value in zip(pairs[:remaining], values[:remaining]):
                    print(
                        json.dumps(
                            {
                                **location(
                                    args.unit,
                                    block,
                                    label,
                                    int(pair),
                                    dataset.ncfs[block - 1],
                                    args.keyorb,
                                ),
                                "coefficient": float(value),
                            }
                        )
                    )
                    emitted += 1
            return 0
        result["elapsed_seconds"] = time.perf_counter() - start
        payload = (
            json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        )
        if args.json:
            protected = {p.resolve() for p in dataset.paths}
            if other is not None:
                protected.update(p.resolve() for p in other.paths)
            protected.update((args.reference / "rcsf.inp",))
            if args.command == "compare":
                protected.add(args.candidate / "rcsf.inp")
            if args.json.resolve() in {p.resolve() for p in protected}:
                raise MCPError("report path would overwrite an input")
            args.json.write_text(payload, encoding="utf-8")
        print(payload, end="")
        return 1 if result["status"] == "different" else 0
    except (OSError, MCPError) as exc:
        print(f"mcp_compare: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
