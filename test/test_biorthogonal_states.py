"""Exercise biorthogonal transforms on one-electron states with a known radial rotation."""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import struct
import subprocess
import tempfile
from itertools import zip_longest
from pathlib import Path

from mixing_stream import Block, read_stream, write_stream
from native_fixture import descriptors, write_radial_stream
from rmcdhf.compare_rwfn import Orbital, read_radial_stream

H = 0.05
RNT = 2e-6
GRID = [RNT * math.expm1(H * i) for i in range(590)]
CASES = {
    "forward": ((1, 2), (1, 2)),
    "reverse": ((2, 1), (2, 1)),
    "extra_left": ((1, 2), (1,)),
    "extra_right": ((1,), (1, 2)),
    "mismatch": ((1, 2), (2, 1)),
    "unordered": ((2, 1, 3), (2, 1, 3)),
    "singular": ((1, 2), (1, 2)),
}
ISOTOPE = (
    "Atomic number:\n1\nMass number:\n1\nFermi a:\n0.5233875553\nFermi c:\n1\n"
    "Mass:\n1\nSpin:\n0.5\nDipole:\n1\nQuadrupole:\n0\n"
)


def orbital_values(n: int, right: bool, singular: bool, r: float) -> float:
    a = 2 * math.exp(-r)
    b = (2 - r) * math.exp(-r / 2) / math.sqrt(8)
    if right:
        a, b = 0.8 * a + 0.6 * b, -0.6 * a + 0.8 * b
    if n == 3:
        return 2 * (27 - 18 * r + 2 * r * r) * math.exp(-r / 3) / (81 * math.sqrt(3))
    return a if n == 1 or singular else b


def state(
    work: Path,
    name: str,
    order: tuple[int, ...],
    singular: bool,
    *,
    swap_csfs: bool = False,
) -> None:
    csf = work / f"{name}.c"
    csf.write_text(
        "Core subshells:\n\nPeel subshells:\n"
        + " ".join(f"{n:3d}s" for n in order)
        + "\nCSF(s):\n"
        + "".join(
            f"  {n}s ( 1)\n      1/2\n       1/2+\n"
            for n in (order[::-1] if swap_csfs else order)
        )
    )
    descriptors(csf)
    csf.unlink()
    right = name == "right"
    rows = []
    for n in order:
        origin = orbital_values(n, right, singular, 0.0) / math.hypot(1, 0.05)
        p = [
            r * orbital_values(n, right, singular, r) / math.hypot(1, 0.05)
            for r in GRID
        ]
        rows.append((n, -1, 0.5 / n**2, origin, p, [0.05 * v for v in p]))
    write_radial_stream(work / f"{name}.wfn", GRID, rows)
    coeff = [1.0] if len(order) == 1 else [0.8, 0.6, 0.0][: len(order)]
    block = Block(1, 1, len(order), [1], [-0.5], struct.pack(f"<{len(order)}d", *coeff))
    write_stream(
        work / f"{name}.mix", 1, len(order), [block], work / f"{name}_desc.parquet"
    )


def overlap(a: Orbital, b: Orbital) -> float:
    # Trapezoidal rule in the uniform log-grid coordinate. Both endpoint terms vanish.
    size = min(len(a.r), len(b.r))
    assert a.r[:size] == b.r[:size]
    return H * math.fsum(
        (r + RNT) * (pa * pb + qa * qb)
        for r, pa, pb, qa, qb in zip(a.r, a.p, b.p, a.q, b.q)
    )


def expansion(
    basis: dict[tuple[int, int], Orbital],
    order: tuple[int, ...],
    vector: tuple[float, ...],
    component: str,
) -> list[float]:
    parts = [getattr(basis[n, -1], component) for n in order]
    return [
        math.fsum(c * v for c, v in zip(vector, row, strict=True))
        for row in zip_longest(*parts, fillvalue=0.0)
    ]


def check_result(
    work: Path, orders: tuple[tuple[int, ...], tuple[int, ...]]
) -> dict[str, float]:
    outputs = []
    coefficients = []
    state_error = 0.0
    for stem, order in zip(("left", "right"), orders, strict=True):
        original = read_radial_stream(work / f"{stem}.wfn")
        result = read_radial_stream(work / f"{stem}.bw")
        assert set(original) == set(result) == {(n, -1) for n in order}
        outputs.append([result[n, -1] for n in order])
        vectors = []
        for suffix in ("mix", "bm"):
            (block,) = read_stream(
                work / f"{stem}.{suffix}", work / f"{stem}_desc.parquet"
            )
            assert (
                block.two_j,
                block.parity,
                block.csf_count,
                block.roots,
                block.energies,
            ) == (1, 1, len(order), [1], [-0.5])
            vectors.append(struct.unpack(f"<{len(order)}d", block.coefficients))
            assert all(math.isfinite(c) for c in vectors[-1])
        coefficients.append((vectors[0], vectors[1]))
        # A one-electron CSF is the orbital itself, so CI compensation must preserve its spinor.
        for n in order:
            a, b = original[n, -1], result[n, -1]
            assert all(
                math.isfinite(v)
                for orb in (a, b)
                for values in (orb.p, orb.q, orb.r)
                for v in values
            )
            assert all(
                math.isclose(r, s, rel_tol=1e-12, abs_tol=1e-18)
                for r, s in zip(a.r, b.r)
            )
        for component in ("p", "q"):
            before = expansion(original, order, vectors[0], component)
            after = expansion(result, order, vectors[1], component)
            state_error = max(
                state_error,
                max(abs(x - y) for x, y in zip_longest(before, after, fillvalue=0.0)),
            )
    gram_error = max(
        abs(overlap(a, b) - int(i == j))
        for i, a in enumerate(outputs[0])
        for j, b in enumerate(outputs[1])
    )
    # Evaluate the identity operator between the two physical states.
    rotation = ((0.8, -0.6), (0.6, 0.8))
    expected = math.fsum(
        coefficients[0][0][i] * rotation[a - 1][b - 1] * coefficients[1][0][j]
        for i, a in enumerate(orders[0])
        for j, b in enumerate(orders[1])
    )
    transformed = math.fsum(
        coefficients[0][1][i] * overlap(a, b) * coefficients[1][1][j]
        for i, a in enumerate(outputs[0])
        for j, b in enumerate(outputs[1])
    )
    return {
        "overlap_error": gram_error,
        "spinor_error": state_error,
        "matrix_element_error": abs(expected - transformed),
    }


def execute(
    work: Path,
    binary: Path,
    launcher: list[str],
    answers: str = "y\nn\nleft\nright\ny\n",
    tag: str = "output",
) -> tuple[int, str]:
    run = subprocess.run(
        [*launcher, str(binary)],
        check=False,
        cwd=work,
        input=answers,
        text=True,
        capture_output=True,
        timeout=45,
        env=os.environ
        | {"MPI_TMP": str(work), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
    )
    log = run.stdout + run.stderr
    (work / f"{tag}.log").write_text(log)
    (work / f"{tag}.exit").write_text(f"{run.returncode}\n")
    return run.returncode, log


def run_case(root: Path, label: str, binary: Path, launcher: list[str]) -> None:
    work = root / label
    work.mkdir()
    (work / "isodata").write_text(ISOTOPE)
    orders = CASES[label]
    for stem, order in zip(("left", "right"), orders, strict=True):
        state(work, stem, order, label == "singular" and stem == "right")
    code, log = execute(work, binary, launcher)
    errors = {
        "mismatch": "initial and final orbital order must match",
        "unordered": "orbital order is neither normal nor reversed",
        "singular": "singular biorthogonal overlap matrix",
    }
    if label in errors:
        assert (
            code != 0
            and errors[label] in log.lower()
            and "Execution complete" not in log
        ), (work, log)
        assert not any(work.glob("*.bm")) and not any(work.glob("*.bw")), work
    else:
        assert code == 0 and "Execution complete" in log, (work, log)
        assert (
            "Reverse orbital ordering"
            if label == "reverse"
            else "Normal orbital ordering"
        ) in log
        assert log.count("NWII") == log.count("NWFF") == 1, (work, log)
        metrics = check_result(work, orders)
        print(label, metrics, flush=True)
        (work / "metrics.json").write_text(json.dumps(metrics) + "\n")
        assert all(math.isfinite(x) and x < 1e-10 for x in metrics.values()), (
            work,
            metrics,
        )


def check_caches(root: Path, binary: Path, launcher: list[str]) -> None:
    work = root / "forward"

    def repeat(tag: str, answers: str = "y\nn\nleft\nright\ny\n") -> str:
        code, log = execute(work, binary, launcher, answers, tag)
        assert code == 0 and "Execution complete" in log, (tag, log)
        return log

    expected = {
        p.name: p.read_bytes()
        for pattern in ("*.bw", "*.bm")
        for p in work.glob(pattern)
    }
    ranks = 2 if launcher else 1
    assert repeat("warm").count("Angular data read from file:") == 2 * ranks
    donor = root / "donor"
    donor.mkdir()
    (donor / "isodata").write_text(ISOTOPE)
    for stem in ("left", "right"):
        state(donor, stem, (1, 2), False, swap_csfs=(stem == "left"))
    code, log = execute(donor, binary, launcher)
    assert code == 0 and "Execution complete" in log, log
    caches = sorted(work.rglob("left.TB"))
    assert len(caches) == ranks
    # Replace a single shard with angular data for a different CSF row order.
    cache = caches[-1]
    shutil.copyfile(donor / cache.relative_to(work), cache)
    assert repeat("wrong-basis").count("Rebuilding incompatible") == 1
    assert all((work / name).read_bytes() == data for name, data in expected.items())
    cache.write_bytes(cache.read_bytes()[:-8])
    assert repeat("truncated").count("Rebuilding incompatible") == 1
    assert all((work / name).read_bytes() == data for name, data in expected.items())
    with cache.open("r+b") as stream:
        stream.seek(4)  # First record's CSF count, after its Fortran record marker.
        stream.write(struct.pack("=i", 0))
    assert repeat("wrong-header").count("Rebuilding incompatible") == 1
    if launcher:
        shutil.copyfile(caches[0], caches[1])
        assert repeat("wrong-rank").count("Rebuilding incompatible") == 1
        assert all(
            (work / name).read_bytes() == data for name, data in expected.items()
        )
    # A request for just J=1/2 uses the same block but a different selection key.
    assert (
        repeat("restricted", "y\nn\nleft\nright\nn\n1\n1\n1\n1\n").count(
            "Rebuilding incompatible"
        )
        == 2 * ranks
    )
    assert all((work / name).read_bytes() == data for name, data in expected.items())
    for stem in ("left", "right"):
        shutil.copyfile(work / f"{stem}.mix", work / f"{stem}.cimix")
    assert (
        repeat("ci", "y\ny\nleft\nright\ny\n").count("Rebuilding incompatible")
        == 2 * ranks
    )
    for stem in ("left", "right"):
        assert (work / f"{stem}.cbm").read_bytes() == expected[f"{stem}.bm"]
    for cache in work.rglob("*.TB"):
        cache.unlink()
    repeat("scratch", "n\nn\nn\nn\nn\nleft\nright\ny\n")
    assert not list(work.rglob("*.TB"))
    assert all((work / name).read_bytes() == data for name, data in expected.items())
    for side, selection in (("initial", "41\n"), ("final", "1\n1\n-1\n")):
        code, log = execute(
            work,
            binary,
            launcher,
            "y\nn\nleft\nright\nn\n" + selection,
            "bad-j-" + side,
        )
        assert code != 0 and "Invalid number of selected J symmetries" in log, log


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--mpiexec")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    launcher = [args.mpiexec, "-n", "2"] if args.mpiexec else []
    with tempfile.TemporaryDirectory(prefix="bio-") as tmp:
        work = args.output.resolve() if args.output else Path(tmp)
        work.mkdir(exist_ok=True)
        invalid = work / "invalid-state"
        invalid.mkdir()
        (invalid / "isodata").write_text(ISOTOPE)
        code, log = execute(invalid, args.binary.resolve(), launcher, "y\nn\nleft\n../right\n")
        assert code != 0 and "Invalid state name" in log, log
        assert log.count("Initial state - Name of state:") == 1, log
        assert log.count("Final state - Name of state:") == 1, log
        assert not list(invalid.glob("*.bw")) and not list(invalid.glob("*.bm"))
        for label in CASES:
            run_case(work, label, args.binary.resolve(), launcher)
        check_caches(work, args.binary.resolve(), launcher)


if __name__ == "__main__":
    main()
