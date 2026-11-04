"""Native GRASP fixtures: scalar/parallel signed doubles, labels and T replay.

No atomic solver or external Python package is needed. Synthetic mixing vectors
are numerical test data, not calculated atomic eigenstates.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile


FIXTURES = Path(__file__).resolve().parent


def prepare(work: Path, case: str, identity: bool = False) -> None:
    work.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(FIXTURES / f"{case}.c", work / "test.c")
    text = (work / "test.c").read_text().splitlines()
    groups: list[list[str]] = [[]]
    for line in text[5:]:
        if line.strip() == "*":
            groups.append([])
        else:
            groups[-1].append(line)
    groups = [g for g in groups if g]
    counts = [len(g) // 3 for g in groups]
    assert all(len(g) % 3 == 0 for g in groups)
    electrons = sum(map(int, re.findall(r"\(\s*(\d+)\)", groups[0][0])))
    nw = len(text[3].split())
    nvec = counts if identity else [min(2, n) for n in counts]

    def ints(*values: int) -> bytes:
        return struct.pack(f"={len(values)}i", *values)

    def doubles(values: list[float]) -> bytes:
        return struct.pack(f"={len(values)}d", *values)

    with (work / "test.m").open("wb") as stream:
        def record(data: bytes) -> None:
            stream.write(ints(len(data)) + data + ints(len(data)))

        record(b"G92MIX")
        record(ints(electrons, sum(counts), nw, sum(nvec),
                    sum(n * v for n, v in zip(counts, nvec)), len(groups)))
        for block, (group, n, v) in enumerate(zip(groups, counts, nvec), 1):
            last = group[2].strip().split()[-1]
            twice_j = int(2 * Fraction(last[:-1]))
            parity = 1 if last[-1] == "+" else -1
            record(ints(block, n, v, twice_j + 1, parity))
            record(ints(*range(1, v + 1)))
            record(doubles([-10.0] + [0.01 * i for i in range(v)]))
            vectors: list[float] = []
            for k in range(v):
                if identity:
                    values = [float(j == k) for j in range(n)]
                else:
                    # All configurations survive default EPSNEW; concentrate
                    # norm in the first configuration to exercise early stop.
                    values = [0.01 * (-1.0 if (j + k) % 2 else 1.0) for j in range(n)]
                    values[k] = math.sqrt(1.0 - (n - 1) * 0.0001)
                vectors.extend(values)
            record(doubles(vectors))
    # Point-nucleus metadata for parser tests only.
    (work / "isodata").write_text(
        "Atomic number:\n10\nMass number:\n20\n"
        "Fermi distribution parameter a:\n0\nFermi distribution parameter c:\n0\n"
        "Mass of nucleus (in amu):\n20\nNuclear spin:\n0\n"
        "Nuclear dipole moment:\n0\nNuclear quadrupole moment:\n0\n"
    )


def input_text(mode: str, unique: bool = False) -> str:
    head = f"test\nn\n{'y' if unique else 'n'}\n"
    if mode == "default":
        return head + "y\n"
    if mode == "omitted":
        return head + "n\ny\n25\n0.000000000001\n0.001\ny\ny\n1\n"
    if mode == "selected":
        # Skip the first and later blocks and reverse two selected level IDs.
        head += "n\nn\n2\n2 1\nn\n"
    else:
        head += "n\ny\n"
    tail = {"write": "y\n", "read": "n\ny\n"}.get(mode, "n\nn\n")
    return head + "0\n0.001\n" + tail + "1\n"


def run(executable: Path, work: Path, mode: str, threads: int,
        force_serial: bool = False, unique: bool = False,
        profile: bool = False, post_serial: bool = False) -> str:
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1", JJ2LSJ_PROFILE="1" if profile else "0",
               JJ2LSJ_POST_SERIAL="1" if post_serial else "0", JJ2LSJ_SERIAL="1" if force_serial else "0")
    result = subprocess.run([str(executable)], input=input_text(mode, unique),
                            text=True, capture_output=True, cwd=work, env=env,
                            timeout=30, check=False)
    output = result.stdout + result.stderr
    (work / "run.log").write_text(output)
    assert result.returncode == 0, output
    assert "Execution complete" in output or "Coefficient dump complete" in output, output
    assert "program stop" not in output.lower(), output
    return output


def vectors(raw: bytes) -> list[list[list[float]]]:
    offset = 0
    blocks = []
    while offset < len(raw):
        _, count, length = struct.unpack_from("=3i", raw, offset)
        offset += 12
        block = []
        for _ in range(count):
            block.append(list(struct.unpack_from(f"={length}d", raw, offset)))
            offset += 8 * length
        blocks.append(block)
    return blocks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("coefficients", type=Path)
    parser.add_argument("--openmp", type=int, choices=(0, 1), required=True)
    args = parser.parse_args()
    executable = args.executable.resolve()
    coefficients = args.coefficients.resolve()
    comparisons = 0
    truncated = False
    with tempfile.TemporaryDirectory(prefix="jj2lsj-regression-") as temporary:
        root = Path(temporary)
        for case in ("even", "odd", "f3"):
            for mode in ("default", "full", "omitted", "selected"):
                scalar = root / f"{case}-{mode}-scalar"
                prepare(scalar, case, identity=(mode in ("full", "selected")))
                run(coefficients, scalar, mode, 4, force_serial=True)
                expected = (scalar / "coefficients.raw").read_bytes()
                assert expected
                parsed = vectors(expected)
                assert all(math.isfinite(x) for b in parsed for v in b for x in v)
                if mode == "full":
                    # Identity jj columns give an independent orthonormality
                    # check of the full signed LS-jj transformation.
                    for block in parsed:
                        for i, left in enumerate(block):
                            for j, right in enumerate(block):
                                product = sum(x*y for x, y in zip(left, right))
                                assert abs(product - float(i == j)) < 2e-12, (case, i, j, product)
                if mode == "omitted":
                    truncated |= any(v[-1] == 0.0 and 0.75 < sum(x*x for x in v) < 0.999999
                                     for b in parsed for v in b)
                run(executable, scalar, mode, 4, force_serial=True)
                if mode == "default":
                    assert (scalar / "test.lsj.lbl").read_bytes() == (
                        FIXTURES / f"{case}.expected.lbl"
                    ).read_bytes(), (case, "upstream label reference")
                outputs = {p.name: p.read_bytes() for p in scalar.glob("test.lsj.*")
                           if p.suffix != ".T"}
                for repeat, threads in enumerate((1, 2, 4, 4)):
                    parallel = root / f"{case}-{mode}-{threads}-{repeat}"
                    prepare(parallel, case, identity=(mode in ("full", "selected")))
                    output = run(coefficients, parallel, mode, threads)
                    assert (parallel / "coefficients.raw").read_bytes() == expected, (case, mode, threads)
                    assert ("jj2lsj OpenMP:" in output) == bool(args.openmp and threads > 1)
                    run(executable, parallel, mode, threads)
                    for name, data in outputs.items():
                        assert (parallel / name).read_bytes() == data, (case, mode, threads, name)
                    comparisons += 1
        assert truncated, "Fixture did not exercise early truncation with a zero LS tail"

        # Sequential T files must stay byte-compatible and replayable even
        # when the executable is launched with multiple OpenMP threads.
        for threads in (1, 4):
            work = root / f"transformation-{threads}"
            prepare(work, "even", identity=True)
            output = run(coefficients, work, "write", threads)
            assert "jj2lsj OpenMP:" not in output
            expected = (work / "coefficients.raw").read_bytes()
            transformation = (work / "test.lsj.T").read_bytes()
            run(coefficients, work, "read", threads)
            assert (work / "coefficients.raw").read_bytes() == expected
            assert (work / "test.lsj.T").read_bytes() == transformation
            if threads == 1:
                scalar_transformation = transformation
            else:
                assert transformation == scalar_transformation

        for threads in (1, 4):
            work = root / f"unique-{threads}"
            prepare(work, "even")
            run(executable, work, "default", threads, unique=True, profile=True)
            labels = [(work / f"test.uni.lsj.{suffix}").read_bytes() for suffix in ("lbl", "sum")]
            assert all(labels)
            with (work / "jj2lsj_omp.profile.csv").open() as handle:
                unique_profile = list(csv.DictReader(handle))
            assert {"unique_read", "unique_assign_summary", "unique_write", "unique_total"} <= {
                row["phase"] for row in unique_profile}
            if threads == 1:
                scalar_labels = labels
            else:
                assert labels == scalar_labels
        # Same block with many ASFs: ordered output despite parallel ranking,
        # and profiling must not change any scientific output bytes.
        for single_post in (True, False):
            work = root / f"profile-{single_post}"
            prepare(work, "even", identity=True)
            run(executable, work, "full", 4, profile=True, post_serial=single_post)
            outputs = {p.name:p.read_bytes() for p in work.glob("test.lsj.*")}
            if single_post:
                post_reference = outputs
            else:
                assert outputs == post_reference
            with (work / "jj2lsj_omp.profile.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            assert rows[-1]["phase"] == "total"
            for phase, suffix in [("label_bytes", ".lsj.lbl"), ("csf_bytes", ".lsj.c"), ("vector_bytes", ".lsj.j")]:
                measured = next(int(r["items"]) for r in rows if r["phase"] == phase)
                assert measured == (work / ("test" + suffix)).stat().st_size
            assert all(float(r["seconds"]) >= 0 for r in rows)
            assert {"csf_read_parse", "jcoup", "mix_read", "ls_rank_wall", "label_format", "label_write"} <= {
                r["phase"] for r in rows}
            workers = [int(r["items"]) for r in rows if r["phase"] == "ls_rank_wall"]
            assert max(workers) == (1 if single_post or not args.openmp else 4)
        assert not (root / "even-full-scalar/jj2lsj_omp.profile.csv").exists()
    print(f"{comparisons} signed-double/file comparisons; orthonormality, truncation, T replay and unique labels passed")


if __name__ == "__main__":
    main()
