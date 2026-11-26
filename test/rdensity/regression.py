"""Exercise real GRASP readers/angular algebra with analytic one-electron states.

All generated inputs/outputs live in a disposable directory. The native rCSFs
CLI prepares descriptors offline; legacy CSF text is removed before every run.
MPI mode also checks the repository's multi-level Ni fixture against the serial
executable on the same suite grid.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from native_fixture import descriptors, write_radial_stream
from mixing_stream import Block, write_stream

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "rmcdhf"))
from compare_rwfn import read_radial_stream


def records(path: Path) -> list[bytes]:
    data = path.read_bytes()
    result: list[bytes] = []
    offset = 0
    while offset < len(data):
        size = struct.unpack_from("=i", data, offset)[0]
        end = offset + 4 + size
        assert size >= 0 and end + 4 <= len(data), "bad Fortran record"
        assert struct.unpack_from("=i", data, end)[0] == size
        result.append(data[offset + 4 : end])
        offset = end + 4
    return result


def write_records(path: Path, rows: list[bytes]) -> None:
    with path.open("wb") as stream:
        for row in rows:
            marker = struct.pack("=i", len(row))
            stream.write(marker + row + marker)


def doubles(values: list[float]) -> bytes:
    return struct.pack(f"={len(values)}d", *values)


def fixture(folder: Path, gap: bool, multiblock: bool, *, native: bool = True) -> None:
    """s state = .8|1s> + .6|ns>; optional pure p_3/2 state.

    The two J weights are 2 and 4. The expected s block is therefore
    outer([.8,.6],[.8,.6])/3 and the p block has occupation 2/3.
    """
    folder.mkdir()
    ns = 3 if gap else 2
    orbitals = [(1, -1), (ns, -1)] + ([(2, -2)] if multiblock else [])
    text = f"Core subshells:\n\nPeel subshells:\n  1s   {ns}s"
    text += "   2p" if multiblock else ""
    text += f"\nCSF(s):\n  1s ( 1)\n      1/2\n       1/2+\n  {ns}s ( 1)\n      1/2\n       1/2+\n"
    if multiblock:
        text += " *\n  2p ( 1)\n      3/2\n       3/2-\n"
    (folder / "case.c").write_text(text)
    # Use a finite nucleus, exercising the suite's default full radial grid.
    (folder / "isodata").write_text(
        "Atomic number:\n1.0\nMass number:\n1.0\nFermi a:\n0.5233875553\n"
        "Fermi c:\n1.0\nMass:\n1.0\nSpin:\n0.5\nDipole:\n1.0\nQuadrupole:\n0.0\n"
    )
    count = len(orbitals)
    blocks = [Block(1, 1, 2, [1], [-0.5], doubles([0.8, 0.6]))]
    if multiblock:
        blocks.append(Block(3, -1, 1, [1], [-0.125], doubles([1.0])))
    if native:
        descriptors(folder / "case.c")
        write_stream(folder / "case.mix", 1, count, blocks, folder / "case_desc.parquet")
    else:
        rows = [b"G92MIX", struct.pack("=6i", 1, count, count, len(blocks), count, len(blocks))]
        for index, block in enumerate(blocks, 1):
            rows += [
                struct.pack("=5i", index, block.csf_count, 1, block.two_j + 1, block.parity),
                struct.pack("=i", 1),
                doubles([block.energies[0], 0.0]),
                block.coefficients,
            ]
        write_records(folder / "case.mix", rows)
    shutil.copyfile(folder / "case.mix", folder / "case.cimix")
    grid = [2e-6 * math.expm1(0.05 * i) for i in range(590)]
    radial_orbitals = []
    for n, kappa in orbitals:
        if kappa == -2:
            p = [r * r * math.exp(-r / 2) / math.sqrt(24) for r in grid]
            pz = 1 / math.sqrt(24)
        elif n == 1:
            p = [2 * r * math.exp(-r) for r in grid]
            pz = 2.0
        else:
            p = [r * (2 - r) * math.exp(-r / 2) / math.sqrt(8) for r in grid]
            pz = 1 / math.sqrt(2)
        radial_orbitals.append((n, kappa, 0.5 / n**2, pz, p, [0.0] * len(grid)))
    if native:
        write_radial_stream(folder / "case.wfn", grid, radial_orbitals)
    else:
        rows = [b"G92RWF"]
        for n, kappa, energy, origin, p, q in radial_orbitals:
            rows += [struct.pack("=iidi", n, kappa, energy, len(grid)), doubles([origin] + p + q), doubles(grid)]
        write_records(folder / "case.wfn", rows)


def run(
    command: list[str], folder: Path, stdin: str | None = None, success: bool = True
) -> str:
    env = dict(
        os.environ,
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
    )
    for text_fixture in folder.glob("*.c"):
        descriptors(text_fixture)
        assert text_fixture.with_name(
            f"{text_fixture.stem}_desc.parquet"
        ).is_file()
        assert text_fixture.with_name(f"{text_fixture.stem}_header.toml").is_file()
        text_fixture.unlink()
    assert not list(folder.glob("*.c")), "rdensity must use descriptor-only CSF inputs"
    result = subprocess.run(
        command,
        check=False,
        cwd=folder,
        input=stdin,
        text=True,
        capture_output=True,
        timeout=90,
        env=env,
    )
    output = result.stdout + result.stderr
    (folder / "run.log").write_text(output)
    if success:
        assert result.returncode == 0 and "Execution complete." in output, output
    else:
        assert result.returncode != 0, output
    return output


def read_blocks(output: str, title: str) -> dict[int, list[list[float]]]:
    lines = output.splitlines()
    kappa = ndim = 0
    result: dict[int, list[list[float]]] = {}
    for i, line in enumerate(lines):
        if line.startswith("KAPPA:"):
            kappa = int(line.split()[-1])
        elif line.startswith("NDIM:"):
            ndim = int(line.split()[-1])
        elif line.strip() == title:
            result[kappa] = [
                [float(x) for x in row.split()[1:]]
                for row in lines[i + 1 : i + ndim + 1]
            ]
    return result


def run_interactive_failure(command: list[str], folder: Path, answers: str) -> str:
    """Give each process a finite stdin without relying on MPI EOF forwarding."""
    input_path = folder / "interactive.inp"
    input_path.write_text(answers)
    launcher = (
        "import os, sys\n"
        "with open(sys.argv[1], 'rb') as stream:\n"
        "    os.dup2(stream.fileno(), 0)\n"
        "os.execv(sys.argv[2], sys.argv[2:])\n"
    )
    # Interactive commands end with the executable and have no application args.
    return run(
        [*command[:-1], sys.executable, "-c", launcher, str(input_path), command[-1]],
        folder,
        success=False,
    )


def close(a: float, b: float, tolerance: float = 2e-10) -> None:
    assert math.isfinite(a) and math.isfinite(b)
    assert abs(a - b) <= tolerance * max(1.0, abs(a), abs(b)), (a, b)


def check_analytic(output: str, folder: Path, multiblock: bool) -> None:
    matrices = read_blocks(output, "density matrix:")
    occupations = read_blocks(output, "occupations in output orbital order:")
    scale = 1 / 3 if multiblock else 1.0
    for actual, expected in zip(
        [value for row in matrices[-1] for value in row],
        [0.64, 0.48, 0.48, 0.36],
        strict=True,
    ):
        close(actual, scale * expected)
    close(occupations[-1][0][0], scale)
    close(occupations[-1][1][0], 0.0)
    if multiblock:
        close(matrices[-2][0][0], 2 / 3)
    close(sum(row[0] for block in occupations.values() for row in block), 1.0)
    # Check the binary NOs reconstruct the weighted radial density. This also
    # exercises the cross term and the P/Q rotations independently of MPI.
    orbitals = read_radial_stream(folder / "case.nw")
    by_kappa: dict[int, int] = {}
    grid = max((orbital.r for orbital in orbitals.values()), key=len)
    reconstructed = [0.0] * len(grid)
    for orbital in orbitals.values():
        slot = by_kappa.get(orbital.kappa, 0)
        occupation = occupations[orbital.kappa][slot][0]
        by_kappa[orbital.kappa] = slot + 1
        for index, (p, q) in enumerate(zip(orbital.p, orbital.q, strict=True)):
            reconstructed[index] += occupation * (p * p + q * q)
    densities: list[dict[float, float]] = []
    for line in (folder / "case.d").read_text().splitlines()[2:]:
        parts = line.split()
        if len(parts) == 3 and parts[0].isdigit():
            densities.append({})
        elif len(parts) == 3:
            r, density, _ = map(float, parts)
            densities[-1][r] = density
    assert len(densities) == (2 if multiblock else 1)
    for i, r in enumerate(grid[1:], 1):
        key = float(f"{r:.10E}")
        expected = densities[0].get(key, 0.0)
        if multiblock:
            expected = (expected + 2 * densities[1].get(key, 0.0)) / 3
        close(reconstructed[i], expected, 3e-9)


def compare_outputs(reference: Path, actual: Path, name: str, ci: bool = False) -> None:
    first = read_radial_stream(reference / f"{name}.nw")
    second = read_radial_stream(actual / f"{name}.nw")
    assert list(first) == list(second)
    for key, orbital in first.items():
        other = second[key]
        assert orbital.energy == other.energy and len(orbital.r) == len(other.r)
        close(orbital.origin, other.origin, 2e-8)
        for values, actual_values in ((orbital.p, other.p), (orbital.q, other.q), (orbital.r, other.r)):
            for x, y in zip(values, actual_values, strict=True):
                close(x, y, 2e-8)
    suffix = ".cd" if ci else ".d"

    def density_rows(
        path: Path,
    ) -> tuple[list[str], dict[tuple[int, float], tuple[float, float]]]:
        headers: list[str] = []
        rows: dict[tuple[int, float], tuple[float, float]] = {}
        for line in path.read_text().splitlines()[2:]:
            parts = line.split()
            if len(parts) == 3 and parts[0].isdigit():
                headers.append(line)
            elif len(parts) == 3:
                r, density, rho = map(float, parts)
                rows[len(headers), r] = (density, rho)
        return headers, rows

    ah, a = density_rows(reference / (name + suffix))
    bh, b = density_rows(actual / (name + suffix))
    assert ah == bh
    # Positive-only output can omit a different subnormal tail row after
    # changing summation order. Missing rows represent zero, within tolerance.
    for key in a.keys() | b.keys():
        for v, w in zip(a.get(key, (0.0, 0.0)), b.get(key, (0.0, 0.0)), strict=True):
            close(v, w, 2e-8)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--serial", required=True, type=Path)
    parser.add_argument("--mpi", type=Path)
    parser.add_argument("--mpiexec")
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--ranks", default=2, type=int)
    parser.add_argument("--data", required=True, type=Path)
    args = parser.parse_args()
    serial = [str(args.serial.resolve())]
    parallel = (
        [args.mpiexec, args.numproc_flag, str(args.ranks), str(args.mpi.resolve())]
        if args.mpi
        else serial
    )
    with tempfile.TemporaryDirectory(prefix="rdensity-regression-") as temporary:
        root = Path(temporary)
        for case, gap, multiblock in [("weighted", False, True), ("gap", True, False)]:
            baseline = root / (case + "-serial")
            fixture(baseline, gap, multiblock)
            output = run(serial + ["case", "--nonci"], baseline)
            check_analytic(output, baseline, multiblock)
            target = root / (case + "-candidate")
            fixture(target, gap, multiblock)
            output = run(parallel + ["case", "--ci"], target)
            # Compare CI output to the same data in the non-CI format.
            shutil.copyfile(target / "case.cd", target / "case.d")
            check_analytic(output, target, multiblock)
            compare_outputs(baseline, target, "case")

        # Existing repository data: 1042 CSFs, 5 J blocks, 9 levels.
        baseline = root / "ni-serial"
        baseline.mkdir()
        for name in ("isodata", "test.mix", "test.wfn", "test_desc.parquet", "test_header.toml"):
            shutil.copyfile(args.data / name, baseline / name)
        ni_output = run(serial + ["test", "--nonci"], baseline)
        occupations = read_blocks(ni_output, "occupations in output orbital order:")
        close(sum(row[0] for block in occupations.values() for row in block), 28.0)
        if args.mpi:
            target = root / "ni-mpi"
            target.mkdir()
            for name in ("isodata", "test.mix", "test.wfn", "test_desc.parquet", "test_header.toml"):
                shutil.copyfile(args.data / name, target / name)
            run(parallel + ["test", "--nonci"], target)
            compare_outputs(baseline, target, "test")
        target = root / "interactive"
        fixture(target, False, False)
        output = run(parallel, target, stdin="y\ncase\nn\n1\n")
        assert output.count("Name of state:") == 1, output
        check_analytic((target / "run.log").read_text(), target, False)
        run(parallel + ["--help"], root)
        run(parallel + ["--unknown"], root, success=False)
        for invalid in ("x" * 24, "x" * 257, "two states"):
            output = run(parallel + [invalid], root, success=False)
            assert "1 to 23 characters" in output, output
        run(parallel + ["missing"], root, success=False)
        if not args.mpi:
            target = root / "archived-mixing"
            fixture(target, False, False)
            for native_suffix, old_suffix, flag in ((".mix", ".m", "--nonci"), (".cimix", ".cm", "--ci")):
                native = target / f"case{native_suffix}"
                archived = native.rename(native.with_suffix(old_suffix))
                output = run(serial + ["case", flag], target, success=False)
                assert native.name in output and not native.exists() and archived.stat().st_size > 0
        # Interactive READ failures must report the runtime message and reach
        # the shared failure path without consuming invalid input or hanging MPI.
        output = run_interactive_failure(parallel, root, "y\n")
        assert "Cannot read state name:" in output, output
        assert "Mixing coefficients from a CI calc.?" not in output, output
        for case, answer in [("missing-order", ""), ("invalid-order", "invalid\n")]:
            target = root / case
            fixture(target, False, False)
            output = run_interactive_failure(parallel, target, "y\ncase\nn\n" + answer)
            assert "Cannot read orbital ordering:" in output, output
            assert not (target / "case.d").exists()
            assert not (target / "case.nw").exists()
        target = root / "unsupported-wavefunction"
        fixture(target, False, False)
        wavefunction = target / "case.wfn"
        previous_stream = b"GRWFRAD1" + wavefunction.read_bytes()[8:]
        for payload in (previous_stream, (args.data / "test.w").read_bytes()):
            wavefunction.write_bytes(payload)
            output = run(parallel + ["case", "--nonci"], target, success=False)
            assert "Invalid radial wavefunction header" in output, output
        # All ranks must terminate on malformed mixing data instead of waiting
        # indefinitely for a collective. Subprocess timeout checks this too.
        target = root / "truncated"
        fixture(target, False, False)
        (target / "case.mix").write_bytes((target / "case.mix").read_bytes()[:-8])
        run(parallel + ["case"], target, success=False)
    print(f"rdensity regression passed ({args.ranks if args.mpi else 0} MPI ranks)")


if __name__ == "__main__":
    main()
