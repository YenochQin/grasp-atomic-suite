"""Synthetic MPI check: tracing must preserve the original radial updates."""
from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str], answers: str, cwd: Path, label: str, trace: str = "0", start_trace: str = "0") -> str:
    env = os.environ | {
        "GRASP_P0_TRACE": trace, "GRASP_START_TRACE": start_trace, "MPI_TMP": str(cwd),
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    }
    result = subprocess.run(command, input=answers, text=True, cwd=cwd, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90, check=False)
    (cwd / f"{label}.log").write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{label}: exit {result.returncode}\n{result.stdout[-4000:]}")
    return result.stdout


def exercise(args: argparse.Namespace, root: Path) -> None:
    (root / "isodata").write_text(
        "Atomic number:\n10\nMass number (integer) :\n20\n"
        "Fermi distribution parameter a:\n0.52338755531043146\n"
        "Fermi distribution parameter c:\n3.0\nMass of nucleus (in amu):\n0\n"
        "Nuclear spin:\n0\nNuclear dipole moment:\n0\nNuclear quadrupole moment:\n0\n"
    )
    (root / "rcsf.inp").write_text(
        "Core subshells:\n\nPeel subshells:\n  1s\nCSF(s):\n  1s ( 2)\n\n         0+\n"
    )
    launcher = [str(args.mpiexec), args.numproc_flag, str(args.ranks)]
    run([str(args.tools / "rwfnestimate")], "y\n2\n*\n", root, "estimate")
    run([*launcher, str(args.tools / "rangular_mpi")], "y\n", root, "angular")
    seed = (root / "rwfn.inp").read_bytes()
    # Explicit grid/accuracy, one EOL level, Method 3 and fixed 0.5 damping.
    answers = (
        "n\nn\ny\nn\ny\n2d-7\n0.05\n0\n590\ny\n1.5625d-8\n"
        "1\n1s\n\n2\ny\nn\nn\nn\ny\nn\n-0.5\nn\nn\nn\nn\n1\n"
    )
    reference: dict[str, bytes] = {}
    for mode in ("0", "1", "all"):
        (root / "rwfn.inp").write_bytes(seed)
        for path in root.glob("rmcdhf.p0.rank*.csv"):
            path.unlink()
        for path in root.glob("rmcdhf.start*.rank*.csv"):
            path.unlink()
        output = run([*launcher, str(args.solver)], answers, root, f"trace-{mode}", mode, mode)
        assert "RMCDHF_MEM_DIAG_MPI: Execution complete." in output
        assert "[P0_TRACE]" not in output, output
        assert "[START_TRACE]" not in output, output
        for name in ("rwfn.out", "rmix.out"):
            data = (root / name).read_bytes()
            assert data
            if mode == "0":
                reference[name] = data
            else:
                assert data == reference[name], f"trace {mode} changed {name}"
        paths = sorted(root.glob("rmcdhf.p0.rank*.csv"))
        expected_ranks = [] if mode == "0" else list(range(args.ranks)) if mode == "all" else [0]
        assert {path.name for path in paths} == {f"rmcdhf.p0.rank{rank}.csv" for rank in expected_ranks}
        for path in paths:
            with path.open() as stream:
                rows = list(csv.DictReader(stream))
            stages = {row["stage"] for row in rows}
            assert {
                "initial", "update_begin", "energy_estimate", "solve_iteration_begin", "before_match",
                "matched_solution", "solve_return", "normalized_solution", "before_damp",
                "damp_mixed", "damp_normalized", "after_damp", "orth_begin",
                "orth_normalized", "orth_end", "cycle_end",
                "homogeneous_start", "homogeneous_out", "homogeneous_in",
                "inhomogeneous_start", "inhomogeneous_out", "inhomogeneous_in",
            } <= stages, stages
            assert all(None not in row for row in rows)
            assert rows[0]["cycle"] == "0" and rows[0]["stage"] == "initial"
            assert rows[-1]["stage"] == "cycle_end"
            assert {row["cycle"] for row in rows} == {"0", "1", "2"}
            rank = rows[0]["rank"]
            with (root / f"rmcdhf.start.rank{rank}.csv").open() as stream:
                start_rows = list(csv.DictReader(stream))
            with (root / f"rmcdhf.start_terms.rank{rank}.csv").open() as stream:
                terms = list(csv.DictReader(stream))
            assert start_rows and terms and all(None not in row for row in start_rows + terms)
            assert {row["itype"] for row in start_rows} == {"1", "2"}
            stage_by_type = {"1": "homogeneous_start", "2": "inhomogeneous_start"}
            for row in start_rows:
                if row["stage"] != "final" or row["point"] == "1":
                    continue
                boundary = next(item for item in rows if item["stage"] == stage_by_type[row["itype"]]
                                and all(item[k] == row[k] for k in ("cycle", "update", "orbital"))
                                and item["inner_iteration"] == row["solve_iteration"])
                scale = float(row["r"]) ** float(row["gamma"])
                for component in ("p", "q"):
                    assert math.isclose(float(boundary[component + row["point"]]),
                                        float(row[component + "_scaled"]) * scale,
                                        rel_tol=1e-13, abs_tol=1e-300)
        for kind in ("start", "start_terms"):
            assert {p.name for p in root.glob(f"rmcdhf.{kind}.rank*.csv")} == {
                f"rmcdhf.{kind}.rank{rank}.csv" for rank in expected_ranks
            }
    print(f"{args.ranks} MPI rank(s): trace off/root/all produce byte-identical rwfn.out and rmix.out")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--tools", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path, required=True)
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--ranks", type=int, default=1)
    args = parser.parse_args()
    if not Path("/proc/meminfo").is_file():
        print("SKIP: upstream memory solver requires Linux /proc/meminfo")
        return 77
    for field in ("solver", "tools", "mpiexec"):
        setattr(args, field, getattr(args, field).resolve())
    with tempfile.TemporaryDirectory(prefix="p0-mpi-", dir="/tmp") as temporary:
        exercise(args, Path(temporary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
