"""Check the installed-command contract with a GRASP case and small G92MIX files."""

from pathlib import Path
import shutil
import struct
import subprocess
import sys
from tempfile import TemporaryDirectory


def records(*payloads: bytes) -> bytes:
    return b"".join(struct.pack("=i", len(data)) + data + struct.pack("=i", len(data)) for data in payloads)


def header(level_count: int, block_count: int) -> bytes:
    return records(b"G92MIX", struct.pack("=6i", 1, 120, 4, level_count, 0, block_count))


def block(
    number: int, j_code: int, parity: int, positions: list[int], shift: float, offsets: list[float]
) -> bytes:
    assert len(positions) == len(offsets)
    count = len(positions)
    csfs = max(positions, default=1)
    block_header = struct.pack("=5i", number, csfs, count, j_code, parity)
    if not count:
        return records(block_header)
    return records(
        block_header,
        struct.pack(f"={count}i", *positions),
        struct.pack(f"={count + 1}d", shift, *offsets),
        bytes(8 * csfs * count),
    )


def label(ordinal: int, j: str, parity: str, energy: float, configuration: str) -> str:
    return (
        f"{ordinal:3d} {j:>4s}     {parity}        {energy:16.9f}     100.000%\n"
        f"       {1.0:12.8f}   {1.0:11.8f}   {configuration}\n"
    )


def run(binary: Path, *files: Path, stdin: str = "") -> subprocess.CompletedProcess[str]:
    return subprocess.run([binary, *files], input=stdin, text=True, capture_output=True, check=False, timeout=5)


def rows(result: subprocess.CompletedProcess[str]) -> list[list[str]]:
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    heading = next(i for i, line in enumerate(lines) if "No Pos  J " in line)
    assert "(a.u.)" in lines[heading + 1] and lines[heading + 1].count("(cm^-1)") == 2
    assert lines[heading + 2].startswith("-----") and lines[-1].startswith("-----")
    return [line.split() for line in lines[heading + 3 : -1]]


def check_real_fixture(binary: Path, data: Path, work: Path) -> None:
    real = rows(run(binary, data / "test.m"))
    assert len(real) == 9
    assert real[0] == ["1", "1", "4", "+", "-1519.3952877", "0.00", "0.00", "3s(2).3p(6).3d(8)3F.4s(2)_3F"]
    assert real[-1] == ["9", "2", "0", "+", "-1519.1147235", "61576.71", "35869.38", "3s(2).3p(6).3d(8)1S.4s(2)_1S"]
    unlabelled = work / "unlabelled.m"
    shutil.copyfile(data / "test.m", unlabelled)
    result = run(binary, stdin=f"{unlabelled}\n\n")
    assert "Configuration" not in result.stdout and "File name ?" in result.stderr
    assert rows(result) == [row[:7] for row in real]
    assert rows(run(binary, stdin=str(unlabelled))) == [row[:7] for row in real]


def check_multiple_files(binary: Path, work: Path) -> None:
    mcdhf = work / ("MCDHF with spaces " + "x" * 76 + ".m")
    rci = work / "RCI.cm"
    # RMCDHF omits blocks 1 and 3; RCI emits a header for its empty block.
    mcdhf.write_bytes(header(3, 4) + block(2, 4, -1, [3, 108], -2.0, [1.0, 1.0]) + block(4, 1, 1, [1], -2.0, [0.0]))
    rci.write_bytes(header(1, 2) + block(1, 999, 999, [], 0.0, []) + block(2, 4, -1, [5], -1.0, [0.0]))
    label_header = " Pos   J   Parity      Energy Total      Comp. of ASF\n"
    long_label = "3s(2)." * 12 + "4s_2S"
    mcdhf.with_suffix(".lsj.lbl").write_text(
        label_header
        + label(1, "3/2", "-", -1.0, long_label)
        + label(2, "3/2", "-", -1.0, "other_degenerate_root")
        + "\n"
        + label(1, "0", "+", 2.0, "opposite_sign"),
        encoding="ascii",
    )
    rci.with_suffix(".lsj.lbl").write_text(label_header + label(1, "3/2", "-", -1.0, "rci_root"), encoding="ascii")
    assert rows(run(binary, mcdhf, rci)) == [
        ["1", "1", "0", "+", "-2.0000000", "0.00", "0.00"],
        ["2", "3", "3/2", "-", "-1.0000000", "219474.63", "219474.63", long_label],
        ["3", "108", "3/2", "-", "-1.0000000", "219474.63", ".00", "other_degenerate_root"],
        ["4", "5", "3/2", "-", "-1.0000000", "219474.63", "0.00", "rci_root"],
    ]


def check_invalid_files(binary: Path, work: Path) -> None:
    bad = work / "invalid.cm"
    for payload, error in (
        (b"", "cannot read G92MIX"),
        (records(b"NOTMIX"), "expected G92MIX"),
        (header(1, 2) + block(2, 1, 1, [1], -1.0, [0.0]), "unexpected symmetry block number"),
        (header(1, 1) + records(struct.pack("=5i", 1, 1, 1, 1, 1), struct.pack("=i", 1)), "block energies"),
        (header(1, 1) + block(1, 1, 1, [1], -1.0, [0.0])[:-5], "mixing coefficient record"),
        (header(0, 0), "no energy levels found"),
    ):
        bad.write_bytes(payload)
        result = run(binary, bad)
        assert result.returncode != 0 and error in result.stderr, result
        assert "No Pos  J " not in result.stdout
    result = run(binary, work / "missing.cm")
    assert result.returncode != 0 and "missing.cm" in result.stderr, result


def main() -> None:
    binary, data = (Path(arg).resolve() for arg in sys.argv[1:])
    with TemporaryDirectory(prefix="rlevels-test-") as directory:
        work = Path(directory)
        check_real_fixture(binary, data, work)
        check_multiple_files(binary, work)
        check_invalid_files(binary, work)


if __name__ == "__main__":
    main()
