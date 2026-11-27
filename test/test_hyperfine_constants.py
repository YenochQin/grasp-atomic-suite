"""Test signed RHFS interaction constants (MHz) and F-dependent matrix elements (Hartree)."""

from __future__ import annotations

import argparse
import math
import shutil
import struct
import tempfile
from fractions import Fraction
from pathlib import Path
from typing import NamedTuple

from mixing_stream import Block, read_stream, write_stream
from test_biorthogonal_states import ISOTOPE
from test_lande_factors import DATA, Level, ls_input, match_levels, run
from test_radiative_transitions import KAYSERS, C

HARTREE_MHZ = KAYSERS * 2.99792458e4
S = Level(1, 1, 1)
P1 = Level(1, -1, 1)
P3 = Level(3, -1, 1)


class Report(NamedTuple):
    diagonal: dict[Level, tuple[float, float]]
    pairs: dict[tuple[Level, Level], tuple[float, float]]
    matrix: dict[tuple[Level, Level, Fraction], float]


def level(fields: list[str]) -> Level:
    root, j, parity = fields
    two_j = 2 * Fraction(j)
    assert two_j.denominator == 1 and parity in ("+", "-"), fields
    return Level(int(two_j), 1 if parity == "+" else -1, int(root))


def pair(a: Level, b: Level) -> tuple[Level, Level]:
    # Only one orientation of the equal-J pair constants is printed.
    return (b, a) if a.two_j == b.two_j and a > b else (a, b)


def read_report(work: Path, ci: bool = False, *, archived: bool = False) -> Report:
    result = Report({}, {}, {})
    for offdiag in (False, True):
        ext = ("cihyperfine" if ci else "hyperfine") + (".off-diagonal" if offdiag else "")
        if archived:
            ext = "hoffd" if offdiag else "h"
        text = (work / f"test.{ext}").read_text()
        assert "A (MHz)" in text and "B (MHz)" in text
        constants, _, matrices = text.split("Interaction constants:")[1].partition(
            "Matrix elements:"
        )
        for line in constants.splitlines():
            fields = line.split()
            if fields and fields[0].isdigit():
                assert len(fields) == 8
                a = level(fields[:3])
                start = 6 if offdiag else 3
                values = tuple(float(x.replace("D", "E")) for x in fields[start:])
                assert all(map(math.isfinite, values)), line
                ab = (values[0], values[1])
                if offdiag:
                    key = pair(a, level(fields[3:6]))
                    assert key not in result.pairs, key
                    result.pairs[key] = ab
                else:
                    assert a not in result.diagonal, a
                    result.diagonal[a] = ab
        if matrices:
            assert offdiag and "Matrix element (a.u.)" in matrices
        for line in matrices.splitlines():
            fields = line.split()
            if fields and fields[0].isdigit():
                assert len(fields) == 8, line
                key_f = (level(fields[:3]), level(fields[3:6]), Fraction(fields[6]))
                value = float(fields[7].replace("D", "E"))
                assert key_f not in result.matrix and math.isfinite(value), line
                result.matrix[key_f] = value
    assert result.diagonal
    for key, ab in result.diagonal.items():
        assert result.pairs[key, key] == ab
    return result


def close(actual: float, expected: float, atol: float) -> None:
    assert math.isclose(actual, expected, rel_tol=2e-10, abs_tol=atol), (
        actual,
        expected,
        atol,
    )


def check_f_energies(report: Report, spin: float) -> None:
    for (bra, ket, f), value in report.matrix.items():
        assert bra.parity == ket.parity
        for state in (bra, ket):
            j = state.two_j / 2
            assert abs(j - spin) <= f <= j + spin
        if bra.two_j == ket.two_j and bra != ket:
            close(value, report.matrix[ket, bra, f], 2e-17)
    # Casimir energy, with K = F(F+1) - I(I+1) - J(J+1).
    for state, (a, b) in report.diagonal.items():
        j = state.two_j / 2
        for two_f in range(
            abs(int(2 * spin) - state.two_j), int(2 * spin) + state.two_j + 1, 2
        ):
            f = Fraction(two_f, 2)
            k = float(f * (f + 1)) - spin * (spin + 1) - j * (j + 1)
            mhz = a * k / 2
            if spin >= 1 and j >= 1:
                mhz += (
                    b
                    * (0.75 * k * (k + 1) - spin * (spin + 1) * j * (j + 1))
                    / (2 * spin * (2 * spin - 1) * j * (2 * j - 1))
                )
            # RHFS omits elements below 1e-14 Hartree.
            value = report.matrix.get((state, state, f))
            if value is None:
                assert abs(mhz / HARTREE_MHZ) <= 1e-14
            else:
                close(value, mhz / HARTREE_MHZ, 2e-17)


def one_electron(work: Path, rhfs: Path, launcher: list[str]) -> None:
    ls_input(work)
    # Hydrogenic Pauli radial integrals and single-electron tensor angular factors.
    magnetic = HARTREE_MHZ * 5.48579909070e-4 / (2 * C**2 * 1.007276466879)
    barn_mhz = HARTREE_MHZ * 1e-24 / (0.52917721067e-8) ** 2
    q = (4 * C) ** -2
    expected_a = {
        S: magnetic * 8 / (3 * (1 + 4 * q)),
        P1: magnetic / (9 * (1 + q)),
        P3: magnetic / (45 * (1 + q)),
    }
    b3 = barn_mhz / 60
    a_off = magnetic / (144 * (1 + q))
    b_off = b3 * math.sqrt(3) / 8 * (1 - 5 * q) / (1 + q)
    expected_pairs = {
        (state, state): (value, b3 if state == P3 else 0.0)
        for state, value in expected_a.items()
    }
    expected_pairs[P3, P1] = (a_off, b_off)
    nuclear_model = ISOTOPE.split("Spin:\n")[0]
    for tag, spin, mu, quadrupole in (
        ("both", 1, 1, 1),
        ("dipole", 1, -2, 0),
        ("quadrupole", 1, 0, 3),
        ("half", 0.5, 1, 0),
        ("spinless", 0, 0, 0),
        ("negative", -1, 0, 0),
    ):
        (work / "isodata").write_text(
            f"{nuclear_model}Spin:\n{spin}\nDipole:\n{mu}\nQuadrupole:\n{quadrupole}\n"
        )
        if spin < 0:
            run(
                work,
                rhfs,
                launcher,
                "y\ntest\nn\n",
                tag,
                error="nonnegative nuclear spin",
            )
            continue
        run(work, rhfs, launcher, "y\ntest\nn\n", tag)
        result = read_report(work)
        assert result.diagonal.keys() == expected_a.keys()
        assert set(result.pairs) == {
            pair(a, b) for a in expected_a for b in expected_a if a.two_j >= b.two_j
        }
        for key, (a, b) in result.pairs.items():
            ma, qb = expected_pairs.get(key, (0.0, 0.0))
            close(a, ma * mu / spin if spin else 0.0, 2e-9)
            close(b, qb * quadrupole, 2e-9)
        check_f_energies(result, spin)
        if not spin:
            assert not result.matrix
        if spin == 1:
            # I=1, J=3/2 -> 1/2 recoupling of the two tensors in the input CSF phase.
            for f, m1, m2 in [
                (Fraction(1, 2), math.sqrt(2), -math.sqrt(200 / 3)),
                (Fraction(3, 2), math.sqrt(5), math.sqrt(20 / 3)),
            ]:
                close(
                    result.matrix[P3, P1, f],
                    (m1 * a_off * mu + m2 * b_off * quadrupole) / HARTREE_MHZ,
                    2e-19,
                )


def vectors(blocks: list[Block]) -> dict[Level, tuple[float, ...]]:
    return {
        Level(b.two_j, b.parity, r): v
        for b in blocks
        for r, v in zip(
            b.roots, struct.iter_unpack(f"<{b.csf_count}d", b.coefficients), strict=True
        )
    }


def align(result: Report, reference: list[Block], target: list[Block]) -> Report:
    matches = match_levels(reference, target)
    left, right = vectors(reference), vectors(target)
    origin = {t: s for s, t in matches.items()}
    sign = {
        t: math.copysign(
            1, math.fsum(x * y for x, y in zip(left[s], right[t], strict=True))
        )
        for s, t in matches.items()
    }
    return Report(
        {origin[t]: ab for t, ab in result.diagonal.items()},
        {
            pair(origin[t], origin[u]): (a * sign[t] * sign[u], b * sign[t] * sign[u])
            for (t, u), (a, b) in result.pairs.items()
        },
        {
            (origin[t], origin[u], f): v * sign[t] * sign[u]
            for (t, u, f), v in result.matrix.items()
        },
    )


def compare(result: Report, reference: Report) -> None:
    for actual, expected in [
        (result.diagonal, reference.diagonal),
        (result.pairs, reference.pairs),
    ]:
        assert actual.keys() == expected.keys()
        for key, (a, b) in actual.items():
            close(a, expected[key][0], 2e-7)
            close(b, expected[key][1], 2e-7)
    assert result.matrix.keys() == reference.matrix.keys()
    for key, v in result.matrix.items():
        close(v, reference.matrix[key], 2e-17)


def nickel(work: Path, rhfs: Path, rci: Path, launcher: list[str]) -> None:
    work.mkdir()
    for name in (
        "isodata",
        "test_desc.parquet",
        "test_header.toml",
        "test.mix",
        "test.wfn",
    ):
        shutil.copyfile(DATA / name, work / name)
    golden = read_report(DATA, archived=True)
    run(work, rhfs, launcher, "y\ntest\nn\n", "mcdhf")
    mcdhf = read_report(work)
    compare(mcdhf, golden)
    check_f_energies(mcdhf, 1)
    descriptor = work / "test_desc.parquet"
    reference = read_stream(work / "test.mix", descriptor)
    # Same fixed-orbital Coulomb Hamiltonian, independently diagonalized by RCI.
    answers = (
        "y\ntest\n"
        + ("n\n" if launcher else "")
        + "n\n" * 6
        + "1 2\n1\n1 2 3\n1\n1 2\n"
    )
    run(work, rci, launcher, answers, "rci")
    path = work / "test.cimix"
    ci = read_stream(path, descriptor)
    run(work, rhfs, launcher, "y\ntest\ny\n", "ci")
    report = read_report(work, True)
    compare(align(report, reference, ci), golden)
    check_f_energies(report, 1)
    # Negate just J=2 root 3, then reverse rows in every block. Both input streams remain present.
    for b in ci:
        row = struct.Struct(f"<{b.csf_count}d")
        rows = [
            row.pack(*(-x for x in v)) if (b.two_j, root) == (4, 3) else row.pack(*v)
            for root, v in zip(b.roots, row.iter_unpack(b.coefficients), strict=True)
        ]
        b.coefficients = b"".join(reversed(rows))
        b.energies.reverse()
        b.roots.reverse()
    write_stream(path, 28, 16, ci, descriptor)
    run(work, rhfs, launcher, "y\ntest\ny\n", "phase")
    changed = read_report(work, True)
    k = (Level(4, 1, 3), Level(2, 1, 1), Fraction(1))
    assert report.matrix[k] != 0 and changed.matrix[k] == -report.matrix[k]
    compare(align(changed, reference, ci), golden)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rhfs", type=Path, required=True)
    parser.add_argument("--rci", type=Path, required=True)
    parser.add_argument("--mpiexec")
    args = parser.parse_args()
    launcher = [args.mpiexec, "-n", "2"] if args.mpiexec else []
    with tempfile.TemporaryDirectory(prefix="rhfs-") as tmp:
        work = Path(tmp)
        one_electron(work / "ls", args.rhfs.resolve(), launcher)
        nickel(work / "ni", args.rhfs.resolve(), args.rci.resolve(), launcher)


if __name__ == "__main__":
    main()
