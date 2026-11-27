"""RIS4 electronic factors and FICAL line shifts from native one-electron and Ni I data."""

from __future__ import annotations

import argparse
import math
import shutil
import struct
import tempfile
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

from mixing_stream import Block, write_stream
from native_fixture import descriptors, write_radial_stream
from test_biorthogonal_states import GRID, ISOTOPE
from test_lande_factors import DATA, ls_input, run
from test_radiative_transitions import C, radial

State = tuple[int, str, str]  # Root, J and parity as printed by RIS4.
AU_TO_GHZ_U = 5.48579909070e-4 * 2 * 109737.31568508 * 29.9792458


@dataclass
class Factors:
    energy: float
    mass: dict[str, tuple[float, ...]] = field(default_factory=dict)
    density: float = 0.0
    fitted: tuple[float, ...] = ()
    corrected: tuple[float, ...] = ()


def numbers(fields: list[str]) -> tuple[float, ...]:
    return tuple(float(x.replace("D", "E")) for x in fields)


def close(
    actual: float, expected: float, *, rel: float = 4e-9, atol: float = 1e-12
) -> None:
    assert math.isclose(actual, expected, rel_tol=rel, abs_tol=atol), (actual, expected)


def state(fields: list[str]) -> State:
    return int(fields[0]), fields[1], fields[2]


def factors(path: Path) -> dict[State, Factors]:
    levels = {}
    section = ""
    lines = iter(path.read_text().splitlines())
    for line in lines:
        columns = line.split()
        if "Level  J Parity" in line:
            section = line.split("Parity", 1)[1].strip()
        if not columns or not columns[0].isdigit():
            continue
        key = state(columns)
        if section == "Energy":
            assert key not in levels
            levels[key] = Factors(numbers(columns[3:4])[0])
        elif "mass shift" in section:
            au = numbers(columns[3:6])
            converted = next(lines).split()
            unit = " ".join(converted[3:])
            assert unit in ("(MHz u)", "(GHz u)")
            ghz = tuple(
                x / (1000 if unit == "(MHz u)" else 1) for x in numbers(converted[:3])
            )
            for x, y in zip(au, ghz, strict=True):
                close(x * AU_TO_GHZ_U, y)
            close(au[0] + au[1], au[2])
            levels[key].mass[section.split()[0]] = ghz
        elif "corrected" in section:
            levels[key].corrected = numbers(columns[3:5])
        elif "electronic factors" in section:
            levels[key].fitted = numbers(columns[3:7])
        else:
            assert "DENS" in section or "density" in section
            levels[key].density = numbers(columns[3:4])[0]
    return levels


def line_tables(path: Path) -> dict[str, dict[tuple[State, State], tuple[float, ...]]]:
    result = {}
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("Upper level"):
            if "NMS" in line:
                label = "mass"
            elif "F0VED0" in line:
                label = "corrected"
            elif "F0" in line:
                label = "field"
            else:
                label = "shift"
            assert label not in result
            rows = result[label] = {}
        cols = line.split()
        if cols and cols[0].isdigit() and len(cols) > 6 and cols[2] in ("+", "-"):
            pair = state(cols), state(cols[3:])
            assert pair not in rows
            rows[pair] = numbers(cols[7:])
    assert set(result) == {"mass", "field", "corrected", "shift"}
    return result


def verify_lines(
    path: Path,
    pairs: dict[tuple[State, State], tuple[Factors, Factors]],
    z: int,
    unit: str,
    column: int,
    moments: tuple[float, ...] = (),
) -> None:
    tables = line_tables(path)
    scale = {"GHz": 1, "MHz": 1000, "meV": 1 / 241.798926}[unit]
    for table in tables.values():
        assert table.keys() == pairs.keys()
    for key, (hi, lo) in pairs.items():
        nms, sms = (
            hi.mass[name][column] - lo.mass[name][column]
            for name in ("Normal", "Specific")
        )
        if hi.fitted:
            f = tuple(x - y for x, y in zip(hi.fitted, lo.fitted, strict=True))
            c0, c1 = (x - y for x, y in zip(hi.corrected, lo.corrected, strict=True))
        else:
            c0 = (hi.density - lo.density) * z * 4.9210908105e-3
            c1 = 0.0
            f = (c0, 0.0, 0.0, 0.0)
        mass = (nms + sms) * (1 / 60 - 1 / 58)
        fs = (
            sum(x * y for x, y in zip(f, moments, strict=True))
            if moments
            else c0 * 0.2 + c1 * 0.2**2
        )
        expected = {
            "mass": (nms, sms, nms + sms),
            "field": f,
            "corrected": (c0, c1),
            "shift": (mass, fs, mass + fs),
        }
        for section, values in expected.items():
            start = (
                1 if section == "mass" else 0
            )  # omit the independent scaling-law estimate
            for actual, target in zip(
                tables[section][key][start : start + len(values)], values, strict=True
            ):
                # In the Fermi case the independently read nuclear moments have six significant digits.
                tolerance = 2e-5 if moments and section == "shift" else 6e-8
                close(actual, target * scale, rel=tolerance, atol=2e-9 * scale)


def two_s_input(work: Path, order: tuple[int, int]) -> None:
    work.mkdir()
    text = work / "test.c"
    text.write_text(
        "Core subshells:\n\nPeel subshells:\n  1s   2s\nCSF(s):\n"
        + "".join(f"  {n}s ( 1)\n      1/2\n       1/2+\n" for n in order)
    )
    descriptors(text)
    (work / "isodata").write_text(ISOTOPE)
    write_radial_stream(work / "test.wfn", GRID, [radial(n, -1) for n in (1, 2)])
    write_stream(
        work / "test.mix",
        1,
        2,
        [Block(1, 1, 2, [1], [-0.5 / order[0] ** 2], struct.pack("<2d", 1, 0))],
        work / "test_desc.parquet",
    )


def one_electron(work: Path, ris: Path, fical: Path, mpi: list[str]) -> None:
    ls_input(work)
    run(work, ris, mpi, "y\ntest\nn\nn\nn\n", "ls")
    levels = factors(work / "test.i")
    for key, n, density in [
        ((1, "1/2", "+"), 1, 1 / math.pi),
        ((1, "1/2", "-"), 2, 9 / (384 * math.pi * C**2)),
        ((1, "3/2", "-"), 2, 0),
    ]:
        s = levels[key]
        close(s.mass["Normal"][0] / AU_TO_GHZ_U, 0.5 / n**2, rel=6e-5)
        assert s.mass["Specific"] == (0.0, 0.0, 0.0)
        close(s.density, density, rel=2e-5, atol=1e-18)
    states = {}
    for name, order in [("lower", (1, 2)), ("upper", (2, 1))]:
        d = work / name
        two_s_input(d, order)
        if name == "upper":
            for ext in ("IOB", "ITB"):
                shutil.copyfile(work / "lower" / f"test.{ext}", d / f"test.{ext}")
        run(d, ris, mpi, "y\ntest\nn\nn\ny\n", "save")
        if name == "upper":
            assert (d / "save.log").read_text().count("incomplete or incompatible") == 2
        result = (d / "test.i").read_bytes()
        run(d, ris, mpi, "y\ntest\nn\nn\n", "reuse")
        assert (d / "test.i").read_bytes() == result
        assert (d / "reuse.log").read_text().count("angular file available") == 2
        shutil.copyfile(d / "test.i", work / (name + ".i"))
        (states[name],) = factors(d / "test.i").values()
    close(states["upper"].density, 1 / (8 * math.pi), rel=4e-5)
    key = (1, "1/2", "+")
    for upper, lower, masses in [
        ("upper", "lower", "60\n58\n0.2"),
        ("lower", "upper", "60\n58\n0.2"),
        ("upper", "lower", "58\n60\n-0.2"),
    ]:
        run(
            work,
            fical,
            [],
            f"n\n{upper}\n{lower}\nn\nn\nboth\nMHz\nn\nn\n{masses}\n",
            "contact-" + upper,
        )
        verify_lines(
            work / f"{upper}.{lower}.fi",
            {(key, key): (states["upper"], states["lower"])},
            1,
            "MHz",
            0,
        )


def reject_damaged_reports(work: Path, fical: Path) -> None:
    lines = (work / "test.i").read_text().rstrip().splitlines()
    last = lines[-1]
    cases = [
        (lines[:4], "reading energy"),
        (
            lines[: next(i for i, s in enumerate(lines) if "Specific mass shift" in s)],
            "Specific mass shift",
        ),
        (
            lines[: next(i for i, s in enumerate(lines) if "Dens." in s) + 1],
            "reading contact density",
        ),
        (lines[:-1] + [last[:40]], "incomplete corrected field factors"),
        (lines[:-1] + [" 999" + last[4:]], "mismatched corrected level"),
        (
            [s.replace("F4 (GHz/fm^6)", "F4 (MHz/fm^6)") for s in lines],
            "expected F4 (GHz/fm^6)",
        ),
    ]
    for index, (damaged, diagnostic) in enumerate(cases):
        (work / "bad.i").write_text("\n".join(damaged) + "\n")
        run(
            work,
            fical,
            [],
            "n\ntest\nbad\nn\ny\npara\nGHz\ny\n",
            f"bad-report-{index}",
            error=diagnostic,
        )
        assert not (work / "test.bad.fi").exists()


def nickel(work: Path, ris: Path, rci: Path, fical: Path, mpi: list[str]) -> None:
    work.mkdir()
    for name in (
        "isodata",
        "test.wfn",
        "test.mix",
        "test_desc.parquet",
        "test_header.toml",
    ):
        shutil.copyfile(DATA / name, work / name)
    for ci in (False, True):
        if ci:
            answers = (
                "y\ntest\n"
                + ("n\n" if mpi else "")
                + "y\ny\n1e-6\n"
                + "n\n" * 4
                + "1 2\n1\n1 2 3\n1\n1 2\n"
            )
            run(work, rci, mpi, answers, "breit")
            for ext in ("IOB", "ITB"):
                (work / f"test.{ext}").unlink()
        mode = "y" if ci else "n"
        run(work, ris, mpi, f"y\ntest\n{mode}\ny\nn\n", "ris-" + mode)
        path = work / ("test.ci" if ci else "test.i")
        levels = factors(path)
        assert len(levels) == 9
        pairs = {}
        for a, b in combinations(levels, 2):
            if levels[b].energy > levels[a].energy:
                a, b = b, a
            pairs[a, b] = levels[a], levels[b]
        for level in levels.values():
            close(level.fitted[0], 28 * 4.9210908105e-3 * level.density, rel=2e-7)
        example = levels[1, "0", "+"]
        golden = NI_REFERENCE[ci]
        for actual, reference in zip(
            (
                *example.mass["Normal"],
                *example.mass["Specific"],
                *example.fitted,
                *example.corrected,
            ),
            golden,
            strict=True,
        ):
            close(actual, reference, rel=5e-8)
        if not ci:
            reject_damaged_reports(work, fical)
            # Save contact-mode coefficients and read them back in Seltzer mode.
            run(work, ris, mpi, "y\ntest\nn\nn\ny\n", "cache-contact")
            run(work, ris, mpi, "y\ntest\nn\ny\n", "cache-seltzer")
            assert (work / "cache-seltzer.log").read_text().count(
                "angular file available"
            ) == 2
            cached = factors(path)
            for key, cold in levels.items():
                warm = cached[key]
                for component in ("Normal", "Specific"):
                    for x, y in zip(
                        cold.mass[component], warm.mass[component], strict=True
                    ):
                        close(x, y)
                for x, y in zip(
                    (*cold.fitted, *cold.corrected, cold.density),
                    (*warm.fitted, *warm.corrected, warm.density),
                    strict=True,
                ):
                    close(x, y)
        for reverse, soph in [
            (False, False),
            (True, False),
            (False, True),
            (True, True),
        ]:
            unit = "meV" if ci else "GHz"
            isotope = "58\n60\n-0.2" if reverse else "60\n58\n0.2"
            if soph:
                nuclei = ["60 3.80 2.3 0 0 0", "58 3.77 2.3 0 0 0"]
                isotope = "\n".join(reversed(nuclei) if reverse else nuclei)
            run(
                work,
                fical,
                [],
                f"n\ntest\ntest\n{mode}\ny\nboth\n{unit}\ny\n"
                + ("y" if soph else "n")
                + f"\n{isotope}\n",
                "fical-" + mode + str((reverse, soph)),
            )
            output = work / ("test.test.cfi" if ci else "test.test.fi")
            moments = tuple(
                numbers(line.split()[-1:])[0]
                for line in output.read_text().splitlines()
                if line.lstrip().startswith("<r^")
            )
            if soph:
                assert len(moments) == 4
                if reverse:
                    moments = tuple(-x for x in moments)
                close(moments[0], 3.80**2 - 3.77**2, rel=2e-5)
            verify_lines(output, pairs, 28, unit, 2, moments)


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in ("ris", "rci", "fical"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--mpiexec")
    args = parser.parse_args()
    mpi = [args.mpiexec, "-n", "2"] if args.mpiexec else []
    with tempfile.TemporaryDirectory(prefix="ris-") as tmp:
        one_electron(Path(tmp) / "ls", args.ris.resolve(), args.fical.resolve(), mpi)
        nickel(
            Path(tmp) / "ni",
            args.ris.resolve(),
            args.rci.resolve(),
            args.fical.resolve(),
            mpi,
        )


# Check serial, fixed Ni I electronic factors at the pre-change baseline.
NI_REFERENCE = {
    False: (
        5759291.253,
        -315483.0638,
        5443808.189,
        -1318965.795,
        43281.22625,
        -1275684.569,
        2696.493498,
        -0.8957368938,
        0.004758816216,
        -1.621990269e-05,
        2669.057915,
        -0.7694974334,
    ),
    True: (
        5759288.221,
        -315482.3903,
        5443805.83,
        -1318963.738,
        43280.40931,
        -1275683.328,
        2696.493524,
        -0.8957369025,
        0.004758816261,
        -1.621990284e-05,
        2669.05794,
        -0.7694974409,
    ),
}

if __name__ == "__main__":
    main()
