"""Exercise the serial and MPI RCI restart entry points with writer-generated files."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile

from mixing_stream import read_stream
from native_fixture import descriptors
from rdensity.regression import records, write_records
from test_csf_descriptors import Record, Shell, text_record


def compare_roots(cold: Path, resumed: Path) -> None:
    before, = read_stream(cold / "case.cimix", cold / "case_desc.parquet")
    after, = read_stream(resumed / "case.cimix", resumed / "case_desc.parquet")
    assert (before.two_j, before.parity, before.csf_count, before.roots) == (0, 1, 4, [1, 2])
    assert (after.two_j, after.parity, after.csf_count, after.roots) == (0, 1, 4, [1, 2])
    assert all(math.isclose(a, b, rel_tol=0.0, abs_tol=1e-10) for a, b in zip(before.energies, after.energies, strict=True))
    for start in (0, 32):
        initial = struct.unpack_from("<4d", before.coefficients, start)
        continued = struct.unpack_from("<4d", after.coefficients, start)
        phase = math.copysign(1.0, sum(a * b for a, b in zip(initial, continued, strict=True)))
        assert max(abs(a - phase * b) for a, b in zip(initial, continued, strict=True)) <= 1e-8


def check_eigenpairs(work: Path, files: list[Path]) -> tuple[float, float]:
    block, = read_stream(work / "case.cimix", work / "case_desc.parquet")
    order = block.csf_count
    matrix = [[0.0] * order for _ in range(order)]
    seen_rows: set[int] = set()
    for file in files:
        shard = records(work / file)
        _, _, orbitals, _ = struct.unpack("=4i", shard[1])
        matrix_header = 7 + 2 * orbitals
        csfs, _, rank, ranks = struct.unpack("=4i", shard[matrix_header])
        assert csfs == order and ranks == len(files)
        rows = range(rank, order, ranks)
        for row, payload in zip(rows, shard[matrix_header + 1 :], strict=True):
            assert row not in seen_rows
            seen_rows.add(row)
            entries, core_shift = struct.unpack_from("=id", payload)
            assert core_shift == 0.0, "fixture uses the unshifted Coulomb matrix"
            values = struct.unpack_from(f"={entries}d", payload, 12)
            partners = struct.unpack_from(f"={entries}i", payload, 12 + 8 * entries)
            assert partners[-1] == row + 1
            for partner, value in zip(partners, values, strict=True):
                assert 1 <= partner <= row + 1 and math.isfinite(value)
                matrix[row][partner - 1] = matrix[partner - 1][row] = value
    assert seen_rows == set(range(order))

    vectors = [struct.unpack_from(f"<{order}d", block.coefficients, 8 * order * i)
               for i in range(len(block.roots))]
    matrix_norm = math.hypot(*(value for row in matrix for value in row))
    residuals = []
    for vector, energy in zip(vectors, block.energies, strict=True):
        residual = math.hypot(*(math.fsum(a * b for a, b in zip(row, vector, strict=True)) - energy * component
                                for row, component in zip(matrix, vector, strict=True)))
        residuals.append(residual / ((matrix_norm + abs(energy)) * math.hypot(*vector)))
    gram_error = max(
        abs(math.fsum(a * b for a, b in zip(left, right, strict=True)) - (i == j))
        for i, left in enumerate(vectors) for j, right in enumerate(vectors)
    )
    assert max(residuals) <= 1e-10 and gram_error <= 1e-12, (work, residuals, gram_error)
    return max(residuals), gram_error


def exercise(root: Path, estimator: Path, rci: Path, launcher: list[str], ranks: int) -> None:
    data = Path(__file__).resolve().parents[1] / "data"
    cold = root / "cold"
    cold.mkdir()
    shutil.copyfile(data / "isodata", cold / "isodata")
    basis = cold / "case.c"
    basis.write_text(
        "Core subshells:\n\nPeel subshells:\n  1s   2s   3s   4s\nCSF(s):\n"
        + "".join(text_record(Record((Shell(f"{n}s", 2),), 0)) for n in range(1, 5))
    )
    descriptors(basis)
    failures = []

    def run(command: list[str], work: Path, answers: str, diagnostic: tuple[str, ...] = ()) -> str:
        try:
            result = subprocess.run(
                command, cwd=work, input=answers, capture_output=True, text=True, timeout=8 if diagnostic else 45,
                env=os.environ | {"MPI_TMP": str(work), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
            )
        except subprocess.TimeoutExpired as error:
            (work / "stdout.log").write_bytes(error.stdout or b"")
            (work / "stderr.log").write_bytes(error.stderr or b"")
            failures.append(f"{work.name}: timed out instead of reporting {diagnostic}")
            return ""
        (work / "stdout.log").write_text(result.stdout)
        (work / "stderr.log").write_text(result.stderr)
        output = result.stdout + result.stderr
        complete = "Execution complete." in output
        if diagnostic:
            if not (result.returncode != 0 and not complete and all(part in output for part in diagnostic)):
                failures.append(f"{work.name}: exit={result.returncode}, expected {diagnostic}; see {work}")
        elif result.returncode != 0 or not complete:
            raise RuntimeError(f"{work.name}: exit={result.returncode}; see {work}\n{output}")
        return output

    run([str(estimator)], cold, "case\ny\n2\n*\n")
    shutil.copyfile(cold / "case.initial.wfn", cold / "case.wfn")
    command = [*launcher, str(rci)]
    run(command, cold, "y\ncase\n" + ("n\n" if launcher else "") + "n\n" * 6 + "1-2\n")
    files = sorted(path.relative_to(cold) for path in cold.rglob("case.rci*.restart"))
    expected = (
        [f"{rank:03d}/case.rci.rank{rank:03d}.restart" for rank in range(ranks)]
        if launcher else ["case.rci.restart"]
    )
    assert [path.as_posix() for path in files] == expected, files
    restart_answers = "n\ncase\ny\nn\n1-2\n"
    warm = root / "warm"
    shutil.copytree(cold, warm)
    (warm / "case.cimix").unlink()
    output = run(command, warm, restart_answers)
    assert output.count("rows read from restart file") == ranks, output
    compare_roots(cold, warm)
    assert all((cold / file).read_bytes() == (warm / file).read_bytes() for file in files)

    # A complete header followed by a clean EOF is a valid unfinished shard.
    partial = root / "partial"
    shutil.copytree(cold, partial)
    (partial / "case.cimix").unlink()
    path = partial / files[-1]
    rows = records(path)
    _, _, orbital_count, _ = struct.unpack("=4i", rows[1])
    matrix_header = 7 + 2 * orbital_count
    write_records(path, rows[:matrix_header])
    run(command, partial, restart_answers)
    compare_roots(cold, partial)
    checks = [check_eigenpairs(work, files) for work in (cold, warm, partial)]
    print(f"{rci.name}/{ranks}: both roots reproduced; max normalized residual "
          f"{max(value[0] for value in checks):.3e}, max orthogonality error "
          f"{max(value[1] for value in checks):.3e}", flush=True)

    cases = [
        ("missing", "cannot open restart file", None, b""),
        ("signature", "restart signature", 0, b"BADRES"),
        ("nuclear-count", "nuclear model parameter count", 3, struct.pack("=i", 3)),
        ("negative-grid", "radial point count", 4, struct.pack("=i", -1)),
        ("oversize-grid", "radial point count", 4, struct.pack("=i", 591)),
        ("truncated-grid", "nuclear potential", 4, struct.pack("=i", 590)),
        ("truncated-header", "radial grid", 6, b""),
        ("orbital-tail", "radial point count for orbital 1", 7, struct.pack("=i", 591)),
        ("block-header", "matrix block header", matrix_header, b""),
        ("matrix-count", "matrix row entry count", matrix_header + 1, struct.pack("=i", 0)),
        ("matrix-row", "matrix row,", matrix_header + 1, struct.pack("=i", 1)),
    ]
    for failed_rank, file in enumerate(files if not launcher or ranks == 2 else []):
        for label, diagnostic, index, payload in cases:
            bad = root / f"bad{failed_rank}-{label}"
            shutil.copytree(cold, bad)
            path = bad / file
            if index is None:
                path.unlink()
            else:
                rows = records(path)
                if label == "orbital-tail":
                    rows[index] = rows[index][:-4] + payload
                elif label in ("nuclear-count", "negative-grid", "oversize-grid", "matrix-count"):
                    rows[index] = payload + rows[index][4:]
                elif label == "truncated-header":
                    rows = rows[:index]
                else:
                    rows[index] = payload
                write_records(path, rows)
            run(command, bad, restart_answers, (diagnostic, file.name, "Input I/O failure"))
            for other in files:
                if other != file:
                    assert (bad / other).read_bytes() == (cold / other).read_bytes(), other
    if failures:
        raise RuntimeError("\n".join(failures))
    if not launcher or ranks == 2:
        print(f"{rci.name}/{ranks}: malformed records and missing shards rejected", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rwfnestimate", type=Path, required=True)
    parser.add_argument("--rci", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--work-dir", type=Path)
    args = parser.parse_args()
    if args.ranks not in (1, 2, 4) or (args.ranks != 1 and args.mpiexec is None):
        parser.error("Choose 1, 2 or 4 ranks and supply an MPI launcher for MPI runs")
    launcher = [str(args.mpiexec.resolve()), args.numproc_flag, str(args.ranks)] if args.mpiexec else []
    if args.work_dir:
        args.work_dir.mkdir()
        exercise(args.work_dir.resolve(), args.rwfnestimate.resolve(), args.rci.resolve(), launcher, args.ranks)
    else:
        with tempfile.TemporaryDirectory(prefix="rci-res-") as temporary:
            exercise(Path(temporary), args.rwfnestimate.resolve(), args.rci.resolve(), launcher, args.ranks)


if __name__ == "__main__":
    main()
