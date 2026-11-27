"""Check producer-written V3 edge cases through serial and MPI applications."""

import argparse
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from mixing_stream import read_stream
from rdensity.regression import records
from native_fixture import descriptors
from test_csf_descriptors import Record, Shell, text_record

ISODATA = Path(__file__).resolve().parents[1] / "data" / "isodata"


def basis(directory: Path, peel: list[str], csfs: list[Record]) -> None:
    directory.mkdir()
    text = directory / "case.c"
    text.write_text(
        "Core subshells:\n\nPeel subshells:\n"
        + "  ".join(peel)
        + "\nCSF(s):\n"
        + "".join(text_record(row) for row in csfs)
    )
    descriptors(text)
    text.unlink()
    shutil.copyfile(ISODATA, directory / "isodata")


def verify(binaries: Path, root: Path, mpiexec: str | None) -> None:
    def run(
        program: str, directory: Path, answers: str, *, failure: str | None = None
    ) -> str:
        binary = binaries / program
        command = [str(binary)]
        if mpiexec and program != "rwfnestimate":
            command = [mpiexec, "-n", "2", str(binary.with_name(program + "_mpi"))]
        result = subprocess.run(
            command,
            cwd=directory,
            input=answers,
            text=True,
            capture_output=True,
            check=False,
            timeout=45,
            env=os.environ
            | {
                "MPI_TMP": str(directory),
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "OMPI_MCA_rmaps_base_oversubscribe": "1",
            },
        )
        output = result.stdout + result.stderr
        if failure is None:
            assert result.returncode == 0 and "Execution complete." in output, output
        else:
            assert result.returncode != 0 and failure in output, output
        return output

    angular_answers = "case\nn\n" if mpiexec else "case\ny\n"
    for orbital, q, two_j, parity in [("5g", 3, 9, 1), ("6h", 11, 11, -1)]:
        work = root / f"{orbital}-{q}"
        shell = Shell(orbital, q, two_j, 1, printed_seniority=1 if q == 3 else None)
        basis(work, [orbital], [Record((shell,), two_j, parity)])
        # Header-only consumers can use the same physical descriptor basis.
        run("rwfnestimate", work, "case\ny\n2\n*\n")
        shutil.copyfile(work / "case.initial.wfn", work / "case.wfn")
        error = f"open j={two_j}/2 shell with occupation {q} exceeds the spin-angular term table"
        for program, answers in [
            ("rangular", angular_answers),
            ("rci", "y\ncase\n" + ("n\n" if mpiexec else "") + "n\n" * 6 + "1\n"),
            ("gj", "y\ncase\nn\n"),
        ]:
            output = run(program, work, answers, failure=error)
            assert "block 1, CSF 1" in output
        assert not list(work.rglob("case.mcp*"))
        assert not (work / "case.gj").exists()

    def check_mcp(work: Path, label: str) -> None:
        run("rangular", work, angular_answers)
        files = sorted(work.rglob("case.mpi.mcpidx" if mpiexec else "case.mcp.30"))
        assert len(files) == (2 if mpiexec else 1), files
        labels = [path.read_bytes()[200:208] if mpiexec else records(path)[2] for path in files]
        assert all(value == label.encode().ljust(8) for value in labels), files

    # A full high-j peel shell follows the closed-shell path in the angular kernel.
    work = root / "full-6h"
    basis(work, ["6h"], [Record((Shell("6h", 12, 0, 0),), 0)])
    check_mcp(work, "   0+")

    # Couple five h^2 shells, g^2 and s^1 to J=115/2. The missing first
    # digit in the old MCP record made this indistinguishable from J=15/2.
    shells = [
        Shell(f"{n}h", 2, 20, 2, 20 * (n - 5) if n > 6 else None) for n in range(6, 11)
    ]
    shells.append(Shell("11g", 2, 16, 2, 116))
    peel = [shell.label for shell in shells] + ["2s", "3s"]
    rows = [Record((*shells, Shell(s, 1, 1, 1)), 115) for s in ("2s", "3s")]
    cases = [("long-label", peel, rows, 115, 1)]

    # Two independent CSFs per block also exercise the two-rank RCI solver.
    rows = [
        Record((Shell("5g", 2, 16, 2), Shell(label, 2, 16, 2)), 32)
        for label in ("6g", "7g")
    ]
    cases.append(("J16", ["5g", "6g", "7g"], rows, 32, 1))
    for label, other, two_j, parity in [
        ("6h", "7h", 11, -1),
        ("10m", "11m", 19, -1),
        ("11n", "12n", 21, 1),
    ]:
        peel = [label, other]
        rows = [
            Record((Shell(orbital, 1, two_j, 1),), two_j, parity) for orbital in peel
        ]
        cases.append((label, peel, rows, two_j, parity))
    peel = ["11n", "12n"]
    rows = [Record((Shell(orbital, 2, 40, 2),), 40) for orbital in peel]
    cases.append(("n-pair", peel, rows, 40, 1))
    for name, peel, rows, two_j, parity in cases:
        work = root / name
        basis(work, peel, rows)
        j_text = f"{two_j}/2" if two_j % 2 else str(two_j // 2)
        sign = "+" if parity == 1 else "-"
        check_mcp(work, f"{j_text:>4}{sign}")
        run("rwfnestimate", work, "case\ny\n2\n*\n")
        shutil.copyfile(work / "case.initial.wfn", work / "case.wfn")
        transverse = "y\ny\n1e-6\n" if name == "n-pair" else "n\nn\n"
        out = run(
            "rci",
            work,
            "y\ncase\n" + ("n\n" if mpiexec else "") + transverse + "n\n" * 4 + "1\n",
        )
        if name == "n-pair":
            assert re.search(
                r"Computing\s+[1-9]\d*\s+Breit integrals of type 2", out
            ), out
        block, = read_stream(work / "case.cimix", work / "case_desc.parquet")
        assert (block.csf_count, block.roots, block.two_j, block.parity) == (len(rows), [1], two_j, parity)
        assert math.isfinite(block.energies[0])
        run("gj", work, "y\ncase\ny\n")
        level = rf"(?m)^\s*1\s+{re.escape(j_text)}\s+{re.escape(sign)}\s+[-+0-9]"
        for suffix in (".cisum", ".cigj"):
            summary = (work / ("case" + suffix)).read_text()
            assert re.search(level, summary), summary
        if not mpiexec and name in ("10m", "11n", "n-pair"):
            run("jj2lsj_omp", work, "case\ny\nn\ny\n")
            labels = (work / "case.ci.lsj.labels").read_text()
            assert re.search(level, labels), labels
            if name in ("10m", "11n"):
                assert f"{name}_2{name[-1].upper()}" in labels, labels


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("binaries", type=Path)
    parser.add_argument("--mpiexec")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="v3-") as tmp:
        verify(args.binaries.resolve(), Path(tmp), args.mpiexec)
