"""Check the installed-command contract with a GRASP case and small native mixing streams."""

from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mixing_stream import Block, write_stream


def block(j_code: int, parity: int, positions: list[int], shift: float, offsets: list[float]) -> Block:
    csfs = max(positions)
    return Block(j_code - 1, parity, csfs, positions, [shift + offset for offset in offsets],
                 bytes(8 * csfs * len(positions)))


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
    mixing = work / "test.mix"
    shutil.copyfile(data / "test.mix", mixing)
    shutil.copyfile(data / "test.lsj.lbl", mixing.with_suffix(".lsj.lbl"))
    result = run(binary, mixing)
    real = rows(result)
    timestamps = r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3} [+-]\d{4}"
    for label in ("Date and Time:", "Finish Date and Time:"):
        assert len(re.findall(rf"^{label} {timestamps}$", result.stderr, re.MULTILINE)) == 1, result.stderr
    assert len(re.findall(r"^Wall time: \d+ seconds$", result.stderr, re.MULTILINE)) == 1, result.stderr
    assert result.stderr.splitlines().count("RLEVELS: Execution complete.") == 1, result.stderr
    assert len(real) == 9
    assert real[0] == ["1", "1", "4", "+", "-1519.3952877", "0.00", "0.00", "3s(2).3p(6).3d(8)3F.4s(2)_3F"]
    assert real[-1] == ["9", "2", "0", "+", "-1519.1147235", "61576.71", "35869.38", "3s(2).3p(6).3d(8)1S.4s(2)_1S"]
    unlabelled = work / "unlabelled.mix"
    shutil.copyfile(mixing, unlabelled)
    result = run(binary, stdin=f"{unlabelled}\n\n")
    assert "Configuration" not in result.stdout and "Name of state:" in result.stderr
    assert rows(result) == [row[:7] for row in real]
    assert rows(run(binary, stdin=str(unlabelled))) == [row[:7] for row in real]


def check_multiple_files(binary: Path, work: Path) -> None:
    mcdhf = work / ("MCDHF with spaces " + "x" * 76 + ".mix")
    rci = mcdhf.with_suffix(".cimix")
    write_stream(mcdhf, 1, 4, [block(4, -1, [3, 108], -2.0, [1.0, 1.0]), block(1, 1, [1], -2.0, [0.0])])
    write_stream(rci, 1, 4, [block(4, -1, [5], -1.0, [0.0])])
    label_header = " Pos   J   Parity      Energy Total      Comp. of ASF\n"
    long_label = "3s(2)." * 12 + "4s_2S"
    mcdhf.with_suffix(".lsj.lbl").write_text("stale label report", encoding="ascii")
    mcdhf.with_suffix(".mcdhf.lsj.labels").write_text(
        label_header
        + label(1, "3/2", "-", -1.0, long_label)
        + label(2, "3/2", "-", -1.0, "other_degenerate_root")
        + "\n"
        + label(1, "0", "+", 2.0, "opposite_sign"),
        encoding="ascii",
    )
    rci.with_suffix(".ci.lsj.labels").write_text(label_header + label(1, "3/2", "-", -1.0, "rci_root"), encoding="ascii")
    assert rows(run(binary, mcdhf, rci)) == [
        ["1", "1", "0", "+", "-2.0000000", "0.00", "0.00"],
        ["2", "3", "3/2", "-", "-1.0000000", "219474.63", "219474.63", long_label],
        ["3", "108", "3/2", "-", "-1.0000000", "219474.63", ".00", "other_degenerate_root"],
        ["4", "5", "3/2", "-", "-1.0000000", "219474.63", "0.00", "rci_root"],
    ]


def check_invalid_files(binary: Path, data: Path, work: Path) -> None:
    bad = work / "invalid.cimix"
    write_stream(bad, 1, 4, [block(1, 1, [1], -1.0, [0.0])])
    valid = bad.read_bytes()
    invalid_row = bytearray(valid)
    struct.pack_into("<q", invalid_row, struct.calcsize("<8s5q64sq"), 0)
    for payload, error in (
        (b"", "Cannot read mixing stream signature"),
        (b"UNKNOWN!", "Expected CSF--MIX"),
        ((data / "test.m").read_bytes(), "Expected CSF--MIX"),
        (b"GRMIX001" + valid[8:], "Expected CSF--MIX"),
        (b"GRMIXTMP" + valid[8:], "Incomplete mixing stream"),
        (valid[:135], "Truncated mixing block directory or level table"),
        (invalid_row, "Invalid mixing block directory entry"),
        (valid[:-5], "Truncated mixing coefficient matrix"),
    ):
        bad.write_bytes(payload)
        result = run(binary, bad)
        assert result.returncode != 0 and error in result.stderr, result
        assert "No Pos  J " not in result.stdout
        assert "RLEVELS: Execution complete." not in result.stderr
    result = run(binary, work / "missing.cimix")
    assert result.returncode != 0 and "missing.cimix" in result.stderr, result


def check_state_selection(binary: Path, work: Path) -> None:
    state = work / "odd.AS1"
    scf = work / "odd.AS1.mix"
    ci = work / "odd.AS1.cimix"
    write_stream(scf, 1, 4, [block(1, 1, [1], -2.0, [0.0])])
    scf_rows = rows(run(binary, scf))
    assert rows(run(binary, state)) == scf_rows
    write_stream(ci, 1, 4, [block(1, 1, [1], -3.0, [0.0])])
    ci_data = ci.read_bytes()
    ci_rows = rows(run(binary, ci))
    assert scf_rows != ci_rows
    for choice, expected in (("n", scf_rows), ("y", ci_rows)):
        result = run(binary, stdin=f"{state}\n{choice}\n\n")
        assert rows(result) == expected and "Mixing coefficients from a CI calc.?" in result.stderr, result
    ambiguous = run(binary, state)
    assert ambiguous.returncode != 0 and str(scf) in ambiguous.stderr and str(ci) in ambiguous.stderr
    scf.unlink()
    assert rows(run(binary, stdin=str(state))) == ci_rows
    for old_suffix in (".m", ".cm"):
        archived = work / f"archived{old_suffix}"
        archived.write_bytes(ci_data)
        result = run(binary, archived)
        assert result.returncode != 0 and "Use a .mix or .cimix" in result.stderr, result
    result = run(binary, work / "archived")
    assert result.returncode != 0 and "cannot open mixing file" in result.stderr, result


def main() -> None:
    binary, data = (Path(arg).resolve() for arg in sys.argv[1:])
    with TemporaryDirectory(prefix="rlevels-test-") as directory:
        work = Path(directory)
        check_real_fixture(binary, data, work)
        check_multiple_files(binary, work)
        check_invalid_files(binary, data, work)
        check_state_selection(binary, work)


if __name__ == "__main__":
    main()
