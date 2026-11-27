"""One-electron E1/M1 transitions from biorthogonalization through the published line list."""

from __future__ import annotations

import argparse
import math
import shutil
import struct
import tempfile
from pathlib import Path
from typing import NamedTuple

from mixing_stream import Block, write_stream
from native_fixture import descriptors, write_radial_stream
from test_biorthogonal_states import GRID, ISOTOPE, RNT, H
from test_descriptor_pipeline import execute

C = 137.035999139
KAYSERS = 2 * 109737.31568508
AU_RATE = (9.10938356e-28 / 1.054571800e-27) * (
    4.803204673e-10**2 / 1.054571800e-27
) ** 2
# Nonrelativistic length-gauge reduced dipoles for 1s,2s -> 2p_3/2.
D1 = math.sqrt(4 / 3) * 256 / (81 * math.sqrt(6))
D2 = -6.0
WEAK_NORM = math.hypot(-D2, 0.997 * D1)
WEAK = (-D2 / WEAK_NORM, 0.997 * D1 / WEAK_NORM)
# Signed reduced elements in the positive input phase. Direct identity-transform
# control at Suite 5a664b58, 0.375 Ha (E1) and 0.005 Ha (M1), report precision.
REFERENCE = {
    (1, "C"): 1.48984,
    (1, "B"): 1.48985,
    (2, "C"): 4.99214e-5,
    (2, "B"): -5.99971,
    (1, "M"): 1.15470,
}


def radial(
    n: int, kappa: int
) -> tuple[int, int, float, float, list[float], list[float]]:
    p, q = [], []
    for r in GRID:
        if n == 1:
            large = 2 * r * math.exp(-r)
            derivative = -large
            origin = 2.0
        elif kappa == -1:
            common = r * math.exp(-r / 2) / math.sqrt(8)
            large = (2 - r) * common
            derivative = (r / 2 - 2) * common
            origin = 1 / math.sqrt(2)
        else:
            common = r * math.exp(-r / 2) / math.sqrt(24)
            large = r * common
            derivative = (3 - r / 2) * common if kappa == 1 else -large / 2
            origin = 1 / math.sqrt(24)
        p.append(large)
        q.append(derivative / (2 * C))
    norm = math.sqrt(
        H
        * math.fsum(
            (r + RNT) * (a * a + b * b) for r, a, b in zip(GRID, p, q, strict=True)
        )
    )
    return (
        n,
        kappa,
        0.5 / n**2,
        origin / norm,
        [x / norm for x in p],
        [x / norm for x in q],
    )


def fixture(work: Path, case: str) -> None:
    work.mkdir()
    (work / "isodata").write_text(ISOTOPE)
    a, b = WEAK if case == "weak" else (1.0, 0.0)
    energy = -0.3125 if case == "half-gap" else -0.5
    odd_energy = {"half-gap": -0.1275, "degenerate": -0.125}.get(case, -0.13)
    even_rows = [
        f"  {ns}s ( 1)\n      1/2\n       1/2+\n"
        for ns in ((2, 1) if case == "donor" else (1, 2))
    ]
    for stem in ("left", "right"):
        if stem == "left":
            peel = "  1s   2s   2p-"
            rows = "".join(even_rows) + " *\n  2p-( 1)\n      1/2\n       1/2-\n"
            orbitals = [(1, -1), (2, -1), (2, 1)]
            blocks = [
                Block(
                    1, 1, 2, [1, 2], [energy, energy], struct.pack("<4d", a, b, -b, a)
                ),
                Block(1, -1, 1, [1], [odd_energy], struct.pack("<d", 1.0)),
            ]
        else:
            peel = "  2p"
            rows = "  2p ( 1)\n      3/2\n       3/2-\n"
            orbitals = [(2, -2)]
            blocks = [
                Block(
                    3,
                    -1,
                    1,
                    [1],
                    [-0.125],
                    struct.pack("<d", -1.0 if case == "flip" else 1.0),
                )
            ]
        csf = work / f"{stem}.c"
        csf.write_text(f"Core subshells:\n\nPeel subshells:\n{peel}\nCSF(s):\n{rows}")
        descriptors(csf)
        csf.unlink()
        write_radial_stream(
            work / f"{stem}.wfn", GRID, [radial(n, k) for n, k in orbitals]
        )
        write_stream(
            work / f"{stem}.mix",
            1,
            len(orbitals),
            blocks,
            work / f"{stem}_desc.parquet",
        )


class Line(NamedTuple):
    energy: float
    rate: float
    gf: float
    strength: float
    element: float | None


class Level(NamedTuple):
    file: str
    root: int
    j: str
    parity: str


def read_report(path: Path, signed: bool) -> dict[tuple[Level, Level, str, str], Line]:
    result = {}
    multipole = ""
    preceding = ""
    for raw in path.read_text().splitlines():
        cols = raw.split()
        if not cols:
            continue
        if cols[0] in ("Electric", "Magnetic"):
            multipole = cols[0][0] + cols[2].split(")")[0]
            preceding = ""
            continue
        if cols[0] in ("f1", "f2", "f") and cols[1] != "=":
            upper = Level(cols[0], int(cols[1]), cols[2], cols[3])
            lower = Level(cols[4], int(cols[5]), cols[6], cols[7])
            energy = float(cols[8])
            gauge = cols[9]
            assert gauge == ("C" if multipole.startswith("E") else "M"), raw
            values = cols[10:]
        elif cols[0] == "B":
            assert preceding == "C", raw
            gauge, values = "B", cols[1:]
        else:
            continue
        fields = [float(x.replace("D", "E")) for x in values]
        assert len(fields) == (4 if signed else 3) and all(
            map(math.isfinite, [energy, *fields])
        ), raw
        key = (upper, lower, multipole, gauge)
        assert key not in result, raw
        result[key] = Line(energy, *fields[:3], fields[3] if signed else None)
        preceding = gauge
    return result


def read_lines(path: Path, reverse: bool, signed: bool) -> dict[tuple[int, str], Line]:
    first, last = ("f1", "f2") if reverse else ("f2", "f1")
    result = {}
    for (upper, lower, multipole, gauge), line in read_report(path, signed).items():
        assert upper == Level(first, 1, "3/2", "-")
        assert multipole in ("E1", "M1")
        assert lower == Level(
            last, lower.root, "1/2", "+" if multipole == "E1" else "-"
        )
        result[lower.root, gauge] = line
    return result


def close(
    actual: float, expected: float, *, rel: float = 1.5e-5, absolute: float = 0
) -> None:
    assert math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute), (
        actual,
        expected,
    )


def check_lines(lines: dict[tuple[int, str], Line], case: str) -> None:
    wanted = set(REFERENCE)
    if case == "degenerate":
        wanted.remove((1, "M"))
    assert set(lines) == wanted, (case, lines)
    for (root, gauge), line in lines.items():
        gap = (0.005 if gauge == "M" else 0.375) / (2 if case == "half-gap" else 1)
        close(line.energy, gap * KAYSERS, rel=0, absolute=0.0051)
        magnetic_factor = 4 * C * C if gauge == "M" else 1.0
        # Upper 2p_3/2 statistical weight is four for both dipoles.
        close(
            line.rate, AU_RATE * (gap / C) ** 3 * line.strength / (3 * magnetic_factor)
        )
        close(line.gf, 2 * gap * line.strength / (3 * magnetic_factor))
        if line.element is not None:
            close(line.element**2, line.strength)
        if case == "half-gap":
            continue
        expected = REFERENCE[root, gauge]
        if case == "weak" and gauge != "M":
            a, b = WEAK
            u, v = REFERENCE[1, gauge], REFERENCE[2, gauge]
            expected = a * u + b * v if root == 1 else -b * u + a * v
        if case == "flip":
            expected = -expected
        # The six printed significant digits of the two component amplitudes
        # limit a prediction for a canceled amplitude to an absolute error.
        error = 6e-6 if case == "weak" else 7e-6 * abs(expected)
        close(math.sqrt(line.strength), abs(expected), rel=0, absolute=error)
        if line.element is not None:
            close(line.element, expected, rel=0, absolute=error)


def regression(root: Path, binary: Path, launcher: list[str]) -> None:
    signed = binary.name == "rtransition_phase"
    bio = binary.with_name("rbiotransform_mpi" if launcher else "rbiotransform")

    def transition(
        work: Path, reverse: bool = False, restricted: bool = False
    ) -> dict[tuple[int, str], Line]:
        first, last = ("right", "left") if reverse else ("left", "right")
        answers = f"y\nn\n{first}\n{last}\nE1,M1\n" + ("y\n" if restricted else "n\n")
        execute(binary, work, answers, launcher, native=True)
        return read_lines(work / f"{first}.{last}.t", reverse, signed)

    all_lines = {}
    for case in ("strong", "weak", "flip", "half-gap", "degenerate", "donor"):
        work = root / case
        fixture(work, case)
        execute(bio, work, "y\nn\nleft\nright\ny\n", launcher, native=True)
        lines = transition(work)
        if case != "donor":
            check_lines(lines, case)
        all_lines[case] = lines
    strong = all_lines["strong"]
    # The finite-wavelength calculation with Pauli small components approaches
    # the hydrogenic E1 dipole and the M1 fine-structure LS limit.
    close(strong[1, "B"].strength, D1**2, rel=2e-4)
    close(strong[2, "B"].strength, D2**2, rel=2e-4)
    close(strong[1, "M"].strength, 4 / 3, rel=2e-4)
    close(all_lines["weak"][1, "B"].strength, 1.92595e-5, rel=0, absolute=2e-10)
    for key, line in strong.items():
        half = all_lines["half-gap"][key]
        close(half.rate / line.rate, (half.strength / line.strength) / 8)
    work = root / "strong"
    assert transition(work, reverse=True) == strong
    caches = sorted(work.rglob("left.right.-1T"))
    ranks = 2 if launcher else 1
    assert len(caches) == ranks

    def repeat(message: str, count: int, restricted: bool = False) -> None:
        assert transition(work, restricted=restricted) == strong
        output = (work / (binary.name + ".stdout")).read_text()
        assert output.count(message) == count, output

    repeat("Angular data read from file:", 2 * ranks)
    shard = caches[-1]
    shutil.copyfile(root / "donor" / shard.relative_to(work), shard)
    repeat("Rebuilding incompatible", 1)
    shard.write_bytes(shard.read_bytes()[:-8])
    repeat("Rebuilding incompatible", 1)
    if launcher:
        shutil.copyfile(caches[0], shard)
        repeat("Rebuilding incompatible", 1)
    repeat("Rebuilding incompatible", ranks, restricted=True)
    # E3 has the correct parity for s -> p but violates the J triangle.
    execute(binary, work, "y\nn\nleft\nright\nE3\n", launcher, native=True)
    report = work / "left.right.t"
    assert "Electric 2**( 3)-pole" in report.read_text()
    assert not read_lines(report, False, signed)
    for cache in work.rglob("left.right.*T"):
        cache.unlink()
    execute(
        binary,
        work,
        "n\nn\nn\nleft\nright\nn\nE1,M1\nn\nKays\nn\nn\n",
        launcher,
        native=True,
    )
    assert read_lines(report, False, signed) == strong
    assert not list(work.rglob("left.right.*T"))


def nickel(work: Path, binary: Path, launcher: list[str]) -> None:
    work.mkdir()
    data = Path(__file__).resolve().parents[1] / "data"
    for source, target in (
        ("isodata", "isodata"),
        ("test.wfn", "ni.bw"),
        ("test.mix", "ni.bm"),
        ("test_desc.parquet", "ni_desc.parquet"),
        ("test_header.toml", "ni_header.toml"),
    ):
        shutil.copyfile(data / source, work / target)
    signed = binary.name == "rtransition_phase"
    report = work / "ni.ni.t"
    execute(binary, work, "y\nn\nni\nni\nM1,E2\nn\nn\n", launcher, native=True)
    lines = read_report(report, signed)
    assert {g: sum(k[-1] == g for k in lines) for g in ("M", "C", "B")} == {
        "M": 14,
        "C": 25,
        "B": 25,
    }
    upper, lower = Level("f", 1, "3", "+"), Level("f", 1, "4", "+")
    for multipole, gauge, rate, element in (
        ("M1", "M", 6.81481e-2, -2.59645),
        ("E2", "C", 6.06554e-8, 0.871642),
        ("E2", "B", 3.03212e-10, 0.0616278),
    ):
        line = lines[upper, lower, multipole, gauge]
        close(line.energy, 1379.17, rel=0, absolute=0.0051)
        close(line.rate, rate)
        close(line.strength, element**2)
        if signed:
            close(line.element, element)
    # Reuse the filenames and change both angular selections. The reference
    # for this run consists exactly of the unequal-J rows from the first run.
    execute(binary, work, "y\nn\nni\nni\nM1,E2\ny\ny\n", launcher, native=True)
    assert read_report(report, signed) == {
        k: v for k, v in lines.items() if k[0].j != k[1].j
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--mpiexec")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    launcher = [args.mpiexec, "-n", "2"] if args.mpiexec else []
    with tempfile.TemporaryDirectory(prefix="rad-") as tmp:
        root = args.output.resolve() if args.output else Path(tmp)
        root.mkdir(exist_ok=True)
        regression(root, args.binary.resolve(), launcher)
        nickel(root / "ni", args.binary.resolve(), launcher)


if __name__ == "__main__":
    main()
