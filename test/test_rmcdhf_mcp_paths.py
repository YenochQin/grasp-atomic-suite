"""Compare the ordinary MPI EOL MCP paths and exercise rank-local format boundaries."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import struct
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from native_fixture import descriptors
from test_csf_descriptors import Record, Shell, text_record

HEADER = struct.Struct("=8s8i12q64s")
BLOCK = struct.Struct("=8s5q")
GROUP = struct.Struct("=4q")
FOOTER = struct.Struct("=8s3q")
ANSWERS = "case\ny\n1\n1\n1\n1s\n1s\n2\n"
ISODATA = Path(__file__).resolve().parents[1] / "data" / "isodata"


def indexed_blocks(path: Path) -> list[tuple[bytes, list[list[tuple[int, int, bytes]]]]]:
    payload = path.read_bytes()
    header = HEADER.unpack_from(payload)
    block_count, channels = header[12], header[13] + 2
    directory = header[-3] - 1
    start = HEADER.size
    blocks = []
    for _ in range(block_count):
        _, _, _, position, columns, rows = BLOCK.unpack_from(payload, start)
        position -= 1
        start += BLOCK.size
        first, counts, totals = (
            struct.unpack_from(f"={channels}q", payload, start + i * 8 * channels)
            for i in range(3)
        )
        start += 24 * channels
        block_channels = []
        for group_start, group_count, total in zip(first, counts, totals, strict=True):
            groups = []
            for index in range(group_start - 1, group_start - 1 + group_count):
                offset, size, label, count = GROUP.unpack_from(payload, directory + index * GROUP.size)
                assert size == count * 16
                groups.append((label, count, payload[offset - 1 : offset - 1 + size]))
            assert sum(group[1] for group in groups) == total
            block_channels.append(groups)
        blocks.append((payload[position : position + 4 * (columns + rows)], block_channels))
    return blocks


def make_basis(folder: Path, *, reversed_order: bool = False) -> None:
    folder.mkdir()
    def pair(orbital: str) -> Record:
        return Record((Shell(orbital, 2),), 0)

    open_pair = Record((Shell("1s", 1, 1, 1), Shell("2s", 1, 1, 1)), 0)
    even = [pair("1s"), open_pair, pair("2s"), pair("2p-")]
    if reversed_order:
        even.reverse()
    odd = Record((Shell("1s", 1, 1, 1), Shell("2p-", 1, 1, 1)), 2, -1)
    path = folder / "case.c"
    path.write_text(
        "Core subshells:\n\nPeel subshells:\n1s 2s 2p-\nCSF(s):\n"
        + "".join(map(text_record, even)) + " *\n" + text_record(odd),
        encoding="ascii",
    )
    descriptors(path)
    shutil.copyfile(ISODATA, folder / "isodata")


def execute(
    command: list[str], work: Path, answers: str = "", *, name: str, failure: str | None = None
) -> str:
    result = subprocess.run(
        command, cwd=work, input=answers, text=True, capture_output=True, timeout=30,
        env=os.environ | {
            "MPI_TMP": str(work),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "OMPI_MCA_rmaps_base_oversubscribe": "1",
        },
    )
    output = result.stdout + result.stderr
    (work / f"{name}.log").write_text(output)
    if failure:
        assert result.returncode != 0 and failure in output, output
        assert "SIGSEGV" not in output, output
    else:
        assert result.returncode == 0 and "Execution complete." in output, output
    return output


def run_paths(root: Path, binaries: Path, mpiexec: str) -> None:
    launcher = [mpiexec, "-n", "2"]
    template = root / "template"
    make_basis(template)
    execute([str(binaries / "rwfnestimate")], template, "case\ny\n2\n*\n", name="orbitals")
    indexed = shutil.copytree(template, root / "indexed")
    execute([*launcher, str(binaries / "rangular_mpi")], indexed, "case\nn\n", name="angular")
    raw = {}
    for rank in ("000", "001"):
        shard = indexed / rank / "case.mpi.mcpidx"
        assert list((indexed / rank).glob("*.mcpidx")) == [shard]
        raw[rank] = indexed_blocks(shard)

    solver = str(binaries / "rmcdhf_mpi")
    stages = set()
    total_reads = 0
    retained_bits = {}
    retained_outputs = {}
    for name in ("scf-first", "scf-repeat"):
        work = shutil.copytree(indexed, root / name)
        log = execute([*launcher, solver, "--mcp-trace", "--mcp-verify"], work, ANSWERS, name="scf")
        assert re.findall(r"Iteration number\s+(\d+)", log) == ["1", "2"], log
        assert len(re.findall(r"^Level +1 +Energy = ", log, re.M)) == 6, log
        for suffix in (".mix", ".wfn"):
            data = (work / ("case" + suffix)).read_bytes()
            assert retained_outputs.setdefault(suffix, data) == data
        run_reads = 0
        for rank in ("000", "001"):
            data = (work / rank / "case.mcp.bits").read_bytes()
            assert retained_bits.setdefault(rank, data) == data
            stages.update(line.split()[0] for line in data.splitlines())
            with (work / rank / "case.mcp.io.csv").open(newline="") as stream:
                trace = list(csv.DictReader(stream))
            calls = [row for row in trace if row["phase"] == "potential_coefficients"]
            assert {row["epoch"] for row in calls} == {"1", "2"}
            assert {row["orbital"] for row in calls} == {"1", "2"}, "frozen multiplier side must be assembled"
            assert Counter((row["epoch"], row["orbital"]) for row in calls)["1", "1"] > 2
            scans = [row for row in trace if row["source"] == "contraction"]
            assert not any(row["source"] == "potential" for row in trace)
            assert max(Counter((row["epoch"], row["block"], row["channel"]) for row in scans).values()) == 1
            for epoch in ("1", "2"):
                used = {
                    (fields[3], fields[4], fields[5])
                    for line in data.splitlines()
                    if (fields := line.split())[0] == b"group" and fields[1] == epoch.encode()
                }
                hits = sum(int(row["groups_read"]) for row in scans if row["epoch"] == epoch)
                assert hits == len(used), (rank, epoch, hits, used)
            run_reads += sum(int(row["bytes_read"]) for row in scans)
            assert not any(row["source"] == "matrix_skip" for row in trace)
            assert {row["block"] for row in trace if row["source"] == "matrix"} == {"1", "2"}
            for row in trace:
                if row["source"] == "matrix" and row["phase"] == "channel":
                    block, channel = int(row["block"]) - 1, int(row["channel"]) - 31
                    groups = raw[rank][block][1][channel]
                    assert int(row["groups_read"]) == len(groups)
                    assert int(row["contributions"]) == sum(group[1] for group in groups)
                    assert int(row["bytes_read"]) == sum(32 + len(group[2]) for group in groups)
        assert run_reads > 0
        total_reads = run_reads
    assert {b"group", b"DA", b"YA", b"XA", b"H", b"multiplier", b"radial"} <= stages
    print(f"Two MPI ranks, two iterations: potential record bytes = {total_reads}")

    memory_solver = binaries / "rmcdhf_mem_mpi"
    if memory_solver.is_file():
        work = shutil.copytree(indexed, root / "memory-scf")
        execute([*launcher, str(memory_solver)], work, ANSWERS, name="scf")
        assert (work / "case.mix").stat().st_size > 0 and (work / "case.wfn").stat().st_size > 0

    reordered = root / "reordered"
    make_basis(reordered, reversed_order=True)
    for name, expected in [
        ("basis", "basis identifier does not match"),
        ("rank", "incompatible rank layout"),
        ("truncated", "MCP indexed truncated shard"),
        ("missing", "MCP indexed open"),
        ("generation", "shards generated by different runs"),
        ("group", "invalid group length"),
        ("sparse-32-limit", "invalid sparse block dimensions or offset"),
        ("sparse-32-overflow", "invalid sparse block dimensions or offset"),
    ]:
        work = shutil.copytree(indexed, root / name)
        shard = work / "001/case.mpi.mcpidx"
        if name == "basis":
            for suffix in ("_desc.parquet", "_header.toml"):
                shutil.copyfile(reordered / f"case{suffix}", work / f"case{suffix}")
        elif name == "rank":
            shutil.copyfile(work / "000/case.mpi.mcpidx", shard)
        elif name == "truncated":
            with shard.open("r+b") as stream:
                stream.truncate(shard.stat().st_size - 1)
        elif name == "missing":
            shard.unlink()
        elif name == "group":
            with shard.open("r+b") as stream:
                header = HEADER.unpack(stream.read(HEADER.size))
                stream.seek(header[-3] - 1 + 8)
                stream.write(struct.pack("=q", 0))
        elif name.startswith("sparse-32-"):
            count = 2**31 - (1 if name.endswith("limit") else 0)
            with shard.open("r+b") as stream:
                stream.seek(HEADER.size + BLOCK.size - 8)
                stream.write(struct.pack("=q", count))
        else:
            with shard.open("r+b") as stream:
                stream.seek(HEADER.size - 64 - 24)
                generation = struct.unpack("=q", stream.read(8))[0] + 1
                stream.seek(-8, os.SEEK_CUR)
                stream.write(struct.pack("=q", generation))
                stream.seek(-24, os.SEEK_END)
                stream.write(struct.pack("=q", generation))
        execute([*launcher, solver], work, ANSWERS, name="scf", failure=expected)
        if name in {"group", "truncated"} and memory_solver.is_file():
            execute([*launcher, str(memory_solver)], work, ANSWERS, name="memory-scf", failure=expected)


def run_format_probe(root: Path, probe: Path, mpiexec: str) -> None:
    def run(scenario: str) -> str:
        output = subprocess.run(
            [mpiexec, "-n", "2", str(probe), scenario], cwd=root, text=True,
            env=os.environ | {"OMPI_MCA_rmaps_base_oversubscribe": "1"},
            capture_output=True, timeout=30,
        )
        assert output.returncode == 0, output.stdout + output.stderr
        return output.stdout

    assert re.search(r"3000000000\s*/\s*5000050000", run("sparsity-count"))
    run("roundtrip")
    path = root / "probe-0.mcpidx"
    data = path.read_bytes()
    header = list(HEADER.unpack_from(data))
    directory = header[-3] - 1
    first, _, label, _ = GROUP.unpack_from(data, directory)
    _, _, second_label, _ = GROUP.unpack_from(data, directory + GROUP.size)
    large_count = 1 + 2**31 // 16
    last_position = first + large_count * 16
    new_directory = last_position + 16
    header[-3] = new_directory
    with path.open("wb") as stream:
        stream.write(HEADER.pack(*header))
        stream.write(data[HEADER.size : first - 1])
        stream.seek(HEADER.size + BLOCK.size + 32)
        stream.write(struct.pack("=q", large_count + 1))
        stream.seek(last_position - 1)
        stream.write(struct.pack("=iiQ", 1, 1, 0x3FF0000000000001))
        stream.write(GROUP.pack(first, large_count * 16, label, large_count))
        stream.write(GROUP.pack(last_position, 16, second_label, 1))
        stream.write(FOOTER.pack(b"MCPDONE1", header[-4], new_directory, 2))
    run("large-offset")
    print("PASS segmented groups, empty rank and channel, 64-bit sparsity and stream positions")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binaries", type=Path)
    parser.add_argument("--mpiexec", required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is None:
        with tempfile.TemporaryDirectory(prefix="rmcdhf-mcp-") as temporary:
            run_paths(Path(temporary), args.binaries.resolve(), args.mpiexec)
            run_format_probe(Path(temporary), args.probe.resolve(), args.mpiexec)
    else:
        args.output.mkdir()
        run_paths(args.output.resolve(), args.binaries.resolve(), args.mpiexec)
        run_format_probe(args.output.resolve(), args.probe.resolve(), args.mpiexec)
