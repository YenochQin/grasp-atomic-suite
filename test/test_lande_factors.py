"""Validate the ordinary, anomalous and total Landé factors through the GJ programs."""

from __future__ import annotations

import argparse
import math
import os
import shutil
import struct
import subprocess
import tempfile
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from typing import NamedTuple

from mixing_stream import Block, read_stream, write_stream
from native_fixture import descriptors, write_radial_stream
from test_biorthogonal_states import GRID, ISOTOPE
from test_radiative_transitions import C, radial

DATA = Path(__file__).resolve().parents[1] / "data"


class Level(NamedTuple):
    two_j: int
    parity: int
    root: int


class Factors(NamedTuple):
    ordinary: float
    anomalous: float
    total: float


def report(path: Path) -> dict[Level, Factors]:
    result = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or not fields[0].isdigit():
            continue
        assert len(fields) == 6 and fields[2] in {"+", "-"}, line
        two_j = 2 * Fraction(fields[1])
        assert two_j.denominator == 1 and two_j > 0, line
        key = Level(int(two_j), 1 if fields[2] == "+" else -1, int(fields[0]))
        assert key not in result, (path, key)
        values = Factors(*(float(x.replace("D", "E")) for x in fields[3:]))
        assert all(map(math.isfinite, values)), line
        assert abs(values.ordinary + values.anomalous - values.total) < 1.1e-10, line
        result[key] = values
    return result


def compare(actual: dict[Level, Factors], expected: dict[Level, Factors]) -> None:
    assert actual.keys() == expected.keys(), (actual.keys(), expected.keys())
    for key, reference in expected.items():
        for value, target, atol in zip(
            actual[key], reference, (5e-10, 5e-12, 5e-10), strict=True
        ):
            assert abs(value - target) < atol, (key, value, target)


def ls_input(work: Path) -> dict[Level, Factors]:
    work.mkdir()
    (work / "isodata").write_text(ISOTOPE)
    csf = work / "test.c"
    csf.write_text(
        "Core subshells:\n\nPeel subshells:\n  1s   2p-  2p\nCSF(s):\n"
        "  1s ( 1)\n      1/2\n       1/2+\n *\n"
        "  2p-( 1)\n      1/2\n       1/2-\n *\n"
        "  2p ( 1)\n      3/2\n       3/2-\n"
    )
    descriptors(csf)
    csf.unlink()
    orbitals = [(1, -1), (2, 1), (2, -2)]
    write_radial_stream(work / "test.wfn", GRID, [radial(n, k) for n, k in orbitals])
    expected = {}
    blocks = []
    for n, kappa in orbitals:
        two_j = 2 * abs(kappa) - 1
        j = two_j / 2
        l = kappa if kappa > 0 else -kappa - 1
        parity = (-1) ** l
        spin_projection = (j * (j + 1) + 0.75 - l * (l + 1)) / (2 * j * (j + 1))
        ls = 1 + spin_projection
        # For these Coulomb Pauli spinors the unnormalized small-component norm is 1/(2nc)^2.
        ordinary = ls / (1 + 1 / (2 * n * C) ** 2)
        delta = 2 * 0.001160 * spin_projection
        expected[Level(two_j, parity, 1)] = Factors(ordinary, delta, ordinary + delta)
        blocks.append(
            Block(two_j, parity, 1, [1], [-0.5 / n**2], struct.pack("<d", 1.0))
        )
    write_stream(work / "test.mix", 1, 3, blocks, work / "test_desc.parquet")
    return expected


def run(
    work: Path,
    binary: Path,
    launcher: list[str],
    answers: str,
    tag: str,
    *,
    error: str = "",
) -> None:
    p = subprocess.run(
        [*launcher, str(binary)],
        input=answers,
        cwd=work,
        text=True,
        capture_output=True,
        check=False,
        timeout=150,
        env=os.environ
        | {
            "MPI_TMP": str(work),
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "OMPI_MCA_rmaps_base_oversubscribe": "1",
        },
    )
    log = p.stdout + p.stderr
    (work / (tag + ".log")).write_text(log)
    if error:
        assert p.returncode != 0 and error in log and "Execution complete" not in log, (
            tag,
            log,
        )
    else:
        assert p.returncode == 0 and "Execution complete" in log, (tag, log)


def match_levels(mcdhf: list[Block], ci: list[Block]) -> dict[Level, Level]:
    """Match complete blocks by maximum absolute CI-vector overlap for these nine Ni I states."""
    matched = {}
    for a, b in zip(mcdhf, ci, strict=True):
        assert (a.two_j, a.parity, a.csf_count) == (b.two_j, b.parity, b.csf_count)
        rows = struct.Struct(f"<{a.csf_count}d")
        overlaps = [
            [
                abs(math.fsum(x * y for x, y in zip(v, w, strict=True)))
                for w in rows.iter_unpack(b.coefficients)
            ]
            for v in rows.iter_unpack(a.coefficients)
        ]
        assert len(a.roots) == len(b.roots)
        order = max(
            permutations(range(len(b.roots))),
            key=lambda seq: sum(overlaps[i][j] for i, j in enumerate(seq)),
        )
        for i, j in enumerate(order):
            assert overlaps[i][j] > 0.999, (a.two_j, a.roots[i], overlaps[i])
            matched[Level(a.two_j, a.parity, a.roots[i])] = Level(
                b.two_j, b.parity, b.roots[j]
            )
    return matched


# Frozen Coulomb+Breit fixed-orbital RCI results for the states of data/test.mix.
# Values are keyed by the original MCDHF root after matching, not by file row.
BREIT = {
    Level(2, 1, 1): Factors(1.4993441093, 1.1596044536e-3, 1.5005037138),
    Level(4, 1, 1): Factors(0.66775135862, -7.6958066102e-4, 0.66698177796),
    Level(4, 1, 2): Factors(1.0428490274, 1.0055929525e-4, 1.0429495867),
    Level(4, 1, 3): Factors(1.4545305661, 1.0555048859e-3, 1.4555860710),
    Level(6, 1, 1): Factors(1.0828552808, 1.9335037238e-4, 1.0830486311),
    Level(8, 1, 1): Factors(1.2493510307, 5.7958981690e-4, 1.2499306205),
    Level(8, 1, 2): Factors(0.99970574031, 4.7364209466e-7, 0.99970621395),
}


def nickel(work: Path, gj: Path, rci: Path, launcher: list[str]) -> None:
    work.mkdir()
    for name in (
        "isodata",
        "test_desc.parquet",
        "test_header.toml",
        "test.mix",
        "test.wfn",
    ):
        shutil.copyfile(DATA / name, work / name)
    run(work, gj, launcher, "y\ntest\nn\n", "mcdhf")
    compare(report(work / "test.gj"), report(DATA / "test.gj"))
    answers = (
        "y\ntest\n"
        + ("n\n" if launcher else "")
        + "y\ny\n1e-6\n"
        + "n\n" * 4
        + "1 2\n1\n1 2 3\n1\n1 2\n"
    )
    run(work, rci, launcher, answers, "rci")
    path = work / "test.cimix"
    descriptor = work / "test_desc.parquet"
    ci = read_stream(path, descriptor)
    identities = match_levels(read_stream(work / "test.mix", descriptor), ci)
    run(work, gj, launcher, "y\ntest\ny\n", "ci")
    reference = report(work / "test.cigj")
    compare(reference, {identities[k]: v for k, v in BREIT.items()})

    # Reverse the storage order, drop the middle J=2 root and flip the ASF phases.
    # The same roots must retain their g factors with both input files present.
    selected = {}
    for block in ci:
        row = struct.Struct(f"<{block.csf_count}d")
        vectors = list(row.iter_unpack(block.coefficients))
        take = [
            i
            for i in reversed(range(len(vectors)))
            if (block.two_j, block.roots[i]) != (4, 2)
        ]
        block.coefficients = b"".join(row.pack(*(-x for x in vectors[i])) for i in take)
        block.energies = [block.energies[i] for i in take]
        block.roots = [block.roots[i] for i in take]
        for root in block.roots:
            key = Level(block.two_j, block.parity, root)
            if key.two_j:
                selected[key] = reference[key]
    write_stream(path, 28, 16, ci, descriptor)
    run(work, gj, launcher, "y\ntest\ny\n", "reordered")
    compare(report(work / "test.cigj"), selected)


def wrong_inputs(work: Path, gj: Path, launcher: list[str]) -> None:
    """Keep the correct descriptor digest and a structurally consistent native stream."""
    work.mkdir()
    for name in ("isodata", "test_desc.parquet", "test_header.toml", "test.wfn"):
        shutil.copyfile(DATA / name, work / name)
    descriptor = work / "test_desc.parquet"
    errors = {
        "electrons": "Mixing dimensions",
        "orbitals": "Mixing dimensions",
        "csfs": "Mixing dimensions",
        "blocks": "Mixing dimensions",
        "partition": "CSF partition",
        "j": "block J or parity",
        "parity": "block J or parity",
        "duplicate": "Duplicate ASF root",
    }
    for case, error in errors.items():
        blocks = read_stream(DATA / "test.mix", descriptor)
        if case == "j":
            blocks[1].two_j = 4
        if case == "parity":
            blocks[1].parity = -1
        if case == "duplicate":
            blocks[2].roots[1] = blocks[2].roots[0]
        if case in ("csfs", "partition", "blocks"):
            if case == "blocks":
                dropped = blocks.pop()
                blocks[-1].csf_count += dropped.csf_count
            else:
                blocks[-1].csf_count -= 1
                if case == "partition":
                    blocks[-2].csf_count += 1
            for block in blocks[-2:]:
                block.coefficients = bytes(8 * len(block.roots) * block.csf_count)
        write_stream(
            work / "test.cimix",
            27 if case == "electrons" else 28,
            15 if case == "orbitals" else 16,
            blocks,
            descriptor,
        )
        run(work, gj, launcher, "y\ntest\ny\n", case, error=error)
        assert not (work / "test.cigj").exists(), case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gj", type=Path, required=True)
    parser.add_argument("--rci", type=Path, required=True)
    parser.add_argument("--mpiexec")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    launcher = [args.mpiexec, "-n", "2"] if args.mpiexec else []
    gj = args.gj.resolve()
    with tempfile.TemporaryDirectory(prefix="lande-") as tmp:
        root = args.output.resolve() if args.output else Path(tmp)
        root.mkdir(exist_ok=True)
        if not launcher:
            debug = root / "debug"
            ls_input(debug)
            run(debug, gj, launcher, "n\ntest\nn\ny\n\ny\nn\nn\nn\nn\nn\nn\n", "debug")
            assert (debug / "gj.dbg").stat().st_size > 0
            assert not list(debug.glob("hyperfine.d*"))
        ls = root / "ls"
        expected = ls_input(ls)
        run(ls, gj, launcher, "y\ntest\nn\n", "ls")
        compare(report(ls / "test.gj"), expected)
        nickel(root / "ni", gj, args.rci.resolve(), launcher)
        wrong_inputs(root / "invalid", gj, launcher)


if __name__ == "__main__":
    main()
