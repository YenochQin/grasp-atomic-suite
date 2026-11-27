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
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from native_fixture import descriptors
from mixing_stream import Block, write_stream
import re
import shutil
import struct
import subprocess
import tempfile


FIXTURES = Path(__file__).resolve().parent


def prepare(work: Path, case: str, identity: bool = False, unused_orbitals: bool = False, *, state: str = "test") -> None:
    work.mkdir(parents=True, exist_ok=True)
    csf_file = work / f"{state}.c"
    shutil.copyfile(FIXTURES / f"{case}.c", csf_file)
    text = csf_file.read_text().splitlines()
    if unused_orbitals:
        text[3] = "  1s " + text[3] + "   5s"
        csf_file.write_text("\n".join(text) + "\n")
    descriptors(csf_file)
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

    blocks = []
    for group, n, v in zip(groups, counts, nvec, strict=True):
        last = group[2].strip().split()[-1]
        twice_j = int(2 * Fraction(last[:-1]))
        parity = 1 if last[-1] == "+" else -1
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
        blocks.append(Block(twice_j, parity, n, list(range(1, v + 1)), [-10.0 + 0.01 * i for i in range(v)],
                            struct.pack(f"<{len(vectors)}d", *vectors)))
    write_stream(work / f"{state}.mix", electrons, nw, blocks, work / f"{state}_desc.parquet")
    # Point-nucleus metadata for parser tests only.
    (work / "isodata").write_text(
        "Atomic number:\n10\nMass number:\n20\n"
        "Fermi distribution parameter a:\n0\nFermi distribution parameter c:\n0\n"
        "Mass of nucleus (in amu):\n20\nNuclear spin:\n0\n"
        "Nuclear dipole moment:\n0\nNuclear quadrupole moment:\n0\n"
    )


def input_text(mode: str, unique: bool = False, *, ci: bool = False, state: str = "test") -> str:
    head = f"{state}\n{'y' if ci else 'n'}\n{'y' if unique else 'n'}\n"
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
        profile: bool = False, post_serial: bool = False, *, ci: bool = False, state: str = "test") -> str:
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1", JJ2LSJ_PROFILE="1" if profile else "0",
               JJ2LSJ_POST_SERIAL="1" if post_serial else "0", JJ2LSJ_SERIAL="1" if force_serial else "0")
    result = subprocess.run([str(executable)], input=input_text(mode, unique, ci=ci, state=state),
                            text=True, capture_output=True, cwd=work, env=env,
                            timeout=30, check=False)
    output = result.stdout + result.stderr
    (work / "run.log").write_text(output)
    assert result.returncode == 0, (executable, work.name, mode, threads, result.returncode, output)
    assert output.count("Name of state:") == 1, output
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
        # Empty peel orbitals stay in the mixing header, but the descriptor
        # reader compacts them before constructing and indexing LS parents.
        reference = root / "compact"
        padded = root / "unused-peel"
        for work, unused in ((reference, False), (padded, True)):
            prepare(work, "f3", unused_orbitals=unused)
            run(executable, work, "default", 1)
            run(coefficients, work, "default", 1)
        for output in ("test.mcdhf.lsj.labels", "coefficients.raw"):
            assert (padded / output).read_bytes() == (reference / output).read_bytes(), output

        # A same-width block with the wrong parity cannot label this basis.
        mixing = padded / "test.mix"
        raw = bytearray(mixing.read_bytes())
        first_block_parity = struct.calcsize("<8s5q64sq")
        old_parity, = struct.unpack_from("<q", raw, first_block_parity)
        struct.pack_into("<q", raw, first_block_parity, -old_parity)
        mixing.write_bytes(raw)
        result = subprocess.run([executable], cwd=padded, input=input_text("default"), text=True,
                                capture_output=True, timeout=30)
        assert result.returncode != 0 and "does not match descriptor symmetry" in result.stderr

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
                    assert (scalar / "test.mcdhf.lsj.labels").read_bytes() == (
                        FIXTURES / f"{case}.expected.lbl"
                    ).read_bytes(), (case, "upstream label reference")
                outputs = {p.name: p.read_bytes() for p in scalar.glob("test.mcdhf.lsj.*")}
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

        # The sequential cache must be deterministic and replayable with different thread counts.
        for threads in (1, 4):
            work = root / f"transformation-{threads}"
            prepare(work, "even", identity=True)
            output = run(coefficients, work, "write", threads)
            assert "jj2lsj OpenMP:" not in output
            expected = (work / "coefficients.raw").read_bytes()
            transformation = (work / "test.lsj.transform.cache").read_bytes()
            run(coefficients, work, "read", threads)
            assert (work / "coefficients.raw").read_bytes() == expected
            assert (work / "test.lsj.transform.cache").read_bytes() == transformation
            if threads == 1:
                scalar_transformation = transformation
            else:
                assert transformation == scalar_transformation

        cache = work / "test.lsj.transform.cache"
        digest_offset = 4 + len(b"LSJTRN02")  # GNU Fortran sequential record marker and cache signature.
        assert set(transformation[digest_offset:digest_offset+64]) <= set(b"0123456789abcdef")
        first_digit = b"a" if transformation[digest_offset:digest_offset+1] != b"a" else b"b"
        foreign = transformation[:digest_offset] + first_digit + transformation[digest_offset+1:]
        for payload, diagnostic in [
            (foreign, "transformation cache CSF basis differs"),
            (transformation[:-16], "Missing JJ2LSJ transformation cache completion record"),
        ]:
            cache.write_bytes(payload)
            result = subprocess.run([str(coefficients)], cwd=work, input=input_text("read"), text=True,
                                    capture_output=True, timeout=30, check=False)
            assert result.returncode != 0 and diagnostic in result.stderr, result.stdout + result.stderr

        # The same state can retain both sets of optional reports.
        work = root / "mcdhf-and-ci"
        prepare(work, "even", identity=True)
        shutil.copyfile(work / "test.mix", work / "test.cimix")
        run(executable, work, "full", 1, unique=True, profile=True)
        suffixes = ("labels", "csf-list", "expansion", "unique-labels", "unique-labels.summary")
        expected_reports = {suffix: (work / f"test.mcdhf.lsj.{suffix}").read_bytes() for suffix in suffixes}
        assert all(expected_reports.values())
        mcdhf_profile = work / "test.mcdhf.lsj.profile.csv"
        expected_profile = mcdhf_profile.read_bytes()
        run(executable, work, "full", 1, unique=True, profile=True, ci=True)
        ci_profile = work / "test.ci.lsj.profile.csv"
        ci_data = ci_profile.read_bytes()
        for suffix, expected_report in expected_reports.items():
            assert (work / f"test.ci.lsj.{suffix}").read_bytes() == expected_report
            assert (work / f"test.mcdhf.lsj.{suffix}").read_bytes() == expected_report
        for legacy_name in ("test.lsj.c", "test.lsj.j", "test.uni.lsj.lbl", "test.uni.lsj.sum"):
            assert not (work / legacy_name).exists()
        prepare(work, "even", state="other-state")
        run(executable, work, "default", 1, profile=True, state="other-state")
        assert mcdhf_profile.read_bytes() == expected_profile
        assert ci_profile.read_bytes() == ci_data
        assert {path.name for path in work.glob("*.profile.csv")} == {
            "test.mcdhf.lsj.profile.csv", "test.ci.lsj.profile.csv", "other-state.mcdhf.lsj.profile.csv",
        }

        for threads in (1, 4):
            work = root / f"unique-{threads}"
            prepare(work, "even")
            run(executable, work, "default", threads, unique=True, profile=True)
            labels = [(work / f"test.mcdhf.lsj.{suffix}").read_bytes()
                      for suffix in ("unique-labels", "unique-labels.summary")]
            assert all(labels)
            with (work / "test.mcdhf.lsj.profile.csv").open() as handle:
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
            profile_path = work / "test.mcdhf.lsj.profile.csv"
            outputs = {p.name: p.read_bytes() for p in work.glob("test.mcdhf.lsj.*") if p != profile_path}
            if single_post:
                post_reference = outputs
            else:
                assert outputs == post_reference
            with profile_path.open() as handle:
                rows = list(csv.DictReader(handle))
            assert rows[-1]["phase"] == "total"
            for phase, suffix in [("label_bytes", "labels"), ("csf_bytes", "csf-list"), ("vector_bytes", "expansion")]:
                measured = next(int(r["items"]) for r in rows if r["phase"] == phase)
                assert measured == (work / f"test.mcdhf.lsj.{suffix}").stat().st_size
            assert all(float(r["seconds"]) >= 0 for r in rows)
            assert {"csf_read_parse", "jcoup", "mix_read", "ls_rank_wall", "label_format", "label_write"} <= {
                r["phase"] for r in rows}
            workers = [int(r["items"]) for r in rows if r["phase"] == "ls_rank_wall"]
            assert max(workers) == (1 if single_post or not args.openmp else 4)
        assert not list((root / "even-full-scalar").glob("*.profile.csv"))
    print(f"{comparisons} signed-double/file comparisons; orthonormality, truncation, T replay and unique labels passed")


if __name__ == "__main__":
    main()
