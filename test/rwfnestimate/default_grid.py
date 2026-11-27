"""Check radial estimates and saved-orbital input on the default grids."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
from tempfile import TemporaryDirectory

from high_l import field, read_radial_stream, run

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from native_fixture import write_radial_stream
from compare_rwfn import Orbital, nodes


def radial_overlap(left: Orbital, right: Orbital, rnt: float, step: float) -> float:
    # Simpson integration in x = log(1 + r/RNT), with dr/dx = r + RNT.
    y = [
        (p1 * p2 + q1 * q2) * (r + rnt)
        for r, p1, p2, q1, q2 in zip(left.r, left.p, right.p, left.q, right.q)
    ]
    y.extend([0.0] * (2 if len(y) % 2 else 3))
    return step / 3 * (y[0] + 4 * math.fsum(y[1:-1:2]) + 2 * math.fsum(y[2:-1:2]))


def verify_orbitals(path: Path, rnt: float, step: float, grid_points: int) -> dict[tuple[int, int], Orbital]:
    orbitals = read_radial_stream(path)
    assert list(orbitals) == [(1, -1), (2, -1), (2, 1), (2, -2)]
    # Use the same quadrature allowance as the existing Ni/Zn high-l regression.
    tolerance = 100 * step**6
    for orbital in orbitals.values():
        assert 1 < len(orbital.r) <= grid_points
        assert all(
            math.isclose(r, rnt * math.expm1(step * i), rel_tol=1e-12)
            for i, r in enumerate(orbital.r)
        ), (path, orbital.n, orbital.kappa)
        assert orbital.energy > 0 and all(math.isfinite(value) for value in (
            orbital.energy, orbital.origin, *orbital.p, *orbital.q
        ))
        angular_l = orbital.kappa if orbital.kappa > 0 else -orbital.kappa - 1
        assert nodes(orbital) == orbital.n - angular_l - 1
        assert abs(radial_overlap(orbital, orbital, rnt, step) - 1.0) < tolerance
    for left, right in combinations(orbitals.values(), 2):
        if left.kappa == right.kappa:
            assert abs(radial_overlap(left, right, rnt, step)) < tolerance
    return orbitals


def write_inputs(directory: Path, point_nucleus: bool) -> None:
    directory.mkdir()
    isodata = (Path(__file__).parent / "fixtures/Ni.isodata").read_text().splitlines()
    if point_nucleus:
        isodata[3] = "0"
    (directory / "isodata").write_text("\n".join(isodata) + "\n")
    peel = "".join(field(n, symbol) for n, symbol in ((1, "s"), (2, "s"), (2, "p-"), (2, "p")))
    header = ["Core subshells:", "", "Peel subshells:", peel, "CSF(s):"]
    (directory / "case_header.toml").write_text(
        "[header_info]\nheader_lines = " + json.dumps(header) + "\nelectron_count = 2\n"
    )


def check_default_grid(binary: Path, root: Path, point_nucleus: bool) -> None:
    rnt = (math.exp(-65 / 16) if point_nucleus else 2e-6) / 28
    step = 0.0625 if point_nucleus else 0.05
    grid_points = 220 if point_nucleus else 590
    for name, option in (("tf", 2), ("hydrogenic", 3)):
        generated = root / name
        write_inputs(generated, point_nucleus)
        run(binary, generated, f"case\ny\n{option}\n*\n")
        reference = verify_orbitals(generated / "case.initial.wfn", rnt, step, grid_points)

        imported = root / f"{name}-import"
        write_inputs(imported, point_nucleus)
        grid = max((orbital.r for orbital in reference.values()), key=len)
        # Reverse the record order and flip orbital phase in the saved input.
        seed = [
            (orbital.n, orbital.kappa, orbital.energy, -orbital.origin,
             [-p for p in orbital.p], [-q for q in orbital.q])
            for orbital in reversed(list(reference.values()))
        ]
        write_radial_stream(imported / "seed.wfn", grid, seed)
        run(binary, imported, "case\ny\n1\nseed.wfn\n*\n")
        restored = verify_orbitals(imported / "case.initial.wfn", rnt, step, grid_points)
        for key, old in reference.items():
            new = restored[key]
            assert new.energy == old.energy
            phase = math.copysign(1.0, radial_overlap(old, new, rnt, step))
            for before, after in ((old.p, new.p), (old.q, new.q)):
                before_padded = before + (0.0,) * (len(grid) - len(before))
                after_padded = after + (0.0,) * (len(grid) - len(after))
                assert max(abs(a - phase * b) for a, b in zip(before_padded, after_padded, strict=True)) < (
                    100 * step**6 * max(abs(value) for value in old.p)
                ), (name, key)


def check_bad_inputs(binary: Path, root: Path) -> None:
    # Use the finite-nucleus default grid and a stream written by rwfnestimate.
    source = (root / "tf/case.initial.wfn").read_bytes()
    _, count, points = struct.unpack_from("<8sii", source)
    short_orbital = bytearray(source)
    invalid_n = bytearray(source)
    too_many_points = bytearray(source)
    duplicate_radius = bytearray(source)
    insufficient_extent = bytearray(source)
    excessive_extent = bytearray(source)
    struct.pack_into("<i", invalid_n, 16, 0)
    struct.pack_into("<i", too_many_points, 24, points + 1)
    struct.pack_into("<i", short_orbital, 24, 1)
    grid_offset = 16 + 32 * count
    struct.pack_into("<d", duplicate_radius, grid_offset + 8, 0.0)
    struct.pack_into("<i", insufficient_extent, 24, 2)
    struct.pack_into("<d", excessive_extent, grid_offset + 8 * (points - 1), 1e7)
    for label, payload, diagnostic in (
        ("invalid-n", invalid_n, "Invalid orbital metadata"),
        ("out-of-range-points", too_many_points, "Invalid orbital metadata"),
        ("truncated", source[:-1], "Invalid radial wavefunction header"),
        ("one-point", short_orbital, "Too few radial points"),
        ("duplicate-radius", duplicate_radius, "strictly increasing"),
        ("insufficient-extent", insufficient_extent, "covers too few active grid points"),
        ("excessive-extent", excessive_extent, "extends beyond the active grid"),
    ):
        directory = root / label
        write_inputs(directory, point_nucleus=False)
        (directory / "seed.wfn").write_bytes(payload)
        result = subprocess.run(
            [str(binary)], cwd=directory, input="case\ny\n1\nseed.wfn\n*\n",
            text=True, capture_output=True, timeout=30, check=False,
        )
        log = result.stdout + result.stderr
        (directory / "rwfnestimate.log").write_text(log)
        assert result.returncode > 0 and diagnostic in log and "Input file: seed.wfn" in log, log
        assert not (directory / "case.initial.wfn").exists()


def check_coarse_input(binary: Path, root: Path) -> None:
    write_inputs(root, point_nucleus=False)
    # Two distinct source radii span many target points on the default grid.
    write_radial_stream(root / "seed.wfn", [0.0, 0.05], [(1, -1, 1.0, 20.0, [0.0, 1.0], [0.0, 0.01])])
    run(binary, root, "case\ny\n1\nseed.wfn\n*\n3\n*\n")
    imported = read_radial_stream(root / "case.initial.wfn")[(1, -1)]
    assert imported.energy == 1.0 and len(imported.p) > 2
    assert all(math.isclose(q, 0.01 * p, rel_tol=1e-12) for p, q in zip(imported.p, imported.q, strict=True))


def exercise(binary: Path, root: Path) -> None:
    for name, point_nucleus in (("point", True), ("finite", False)):
        directory = root / name
        directory.mkdir()
        check_default_grid(binary, directory, point_nucleus)
        print(f"{name}: verified Thomas–Fermi and hydrogenic orbitals and saved-orbital import")
    check_coarse_input(binary, root / "coarse")
    check_bad_inputs(binary, root / "finite")
    print("Rejected malformed saved-orbital streams before publishing initial orbitals")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--workdir", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if args.workdir:
        args.workdir.mkdir()
        exercise(binary, args.workdir)
    else:
        with TemporaryDirectory(prefix="rwfn-default-") as directory:
            exercise(binary, Path(directory))


if __name__ == "__main__":
    main()
