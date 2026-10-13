#!/usr/bin/env python3
"""Small synthetic two-electron test of strict stopping and .w restart I/O."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from compare_rwfn import read_rwfn


def inputs(cycles: int) -> str:
    # Explicit finite-nucleus grid; one EOL level has no weight-menu prompt.
    return (
        "n\nn\ny\nn\ny\n2d-7\n0.05\n0\n590\ny\n1.5625d-8\n"
        f"1\n1s\n\n{cycles}\nn\n1\n"
    )


def run(command: list[str], stdin: str, cwd: Path, label: str) -> str:
    env = os.environ.copy()
    env.update(MPI_TMP=str(cwd), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    result = subprocess.run(
        command, input=stdin, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, cwd=cwd, env=env, timeout=90, check=False,
    )
    (cwd / f"{label}.log").write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{label}: exit {result.returncode}\n{result.stdout[-4000:]}")
    return result.stdout


def compare_wavefunctions(left: Path, right: Path) -> None:
    old, new = read_rwfn(left), read_rwfn(right)
    assert old.keys() == new.keys() == {(1, -1)}
    a, b = old[1, -1], new[1, -1]
    common = min(len(a.r), len(b.r))
    assert a.r[:common] == b.r[:common]
    assert abs(a.energy - b.energy) < 1e-10
    assert all(math.isfinite(x) for x in (*b.p, *b.q))
    for left_component, right_component in ((a.p, b.p), (a.q, b.q)):
        size = max(len(left_component), len(right_component))
        padded_left = left_component + (0.0,) * (size - len(left_component))
        padded_right = right_component + (0.0,) * (size - len(right_component))
        assert max(abs(x-y) for x, y in zip(padded_left, padded_right)) < 1e-9


def exercise(args: argparse.Namespace, root: Path) -> None:
    # Synthetic Z=10 finite nucleus with 1s^2. No external research data.
    (root / "isodata").write_text(
        "Atomic number:\n10\nMass number (integer) :\n20\n"
        "Fermi distribution parameter a:\n0.52338755531043146\n"
        "Fermi distribution parameter c:\n3.0\n"
        "Mass of nucleus (in amu):\n0\nNuclear spin:\n0\n"
        "Nuclear dipole moment:\n0\nNuclear quadrupole moment:\n0\n"
    )
    (root / "rcsf.inp").write_text(
        "Core subshells:\n\nPeel subshells:\n  1s\nCSF(s):\n  1s ( 2)\n\n         0+\n"
    )
    launcher = [str(args.mpiexec), args.numproc_flag, str(args.ranks)]
    run([str(args.tools / "rwfnestimate")], "y\n2\n*\n", root, "estimate")
    run([*launcher, str(args.tools / "rangular_mpi")], "y\n", root, "angular")
    strict = [*launcher, str(args.strict)]
    output = run(strict, inputs(1), root, "limited")
    assert "SCF termination: ITERATION_LIMIT (not converged)" in output
    with (root / "rwfn.out.scf.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1 and rows[0]["converged"] == "F"
    assert rows[0]["energy_ok"] == "F"
    for row in rows:
        assert math.isfinite(float(row["weighted_energy"]))
    shutil.copy2(root / "rwfn.out", root / "restart.w")
    assert (root / "rwfn.out.iter000000").is_file()
    compare_wavefunctions(root / "restart.w", root / "rwfn.out.iter000001")

    # Same first sweep as the existing suite memory solver: stopping and
    # diagnostics must not perturb the radial update itself.
    run([*launcher, str(args.baseline)], inputs(1), root, "baseline")
    compare_wavefunctions(root / "restart.w", root / "rwfn.out")

    shutil.copy2(root / "restart.w", root / "rwfn.inp")
    output = run(strict, inputs(100), root, "restart")
    assert "SCF termination: STRICT_CONVERGED" in output
    compare_wavefunctions(root / "restart.w", root / "rwfn.out.iter000000")
    with (root / "rwfn.out.scf.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert 3 <= len(rows) <= 100
    assert rows[0]["energy_ok"] == "F" and rows[0]["pass_streak"] == "0"
    assert [i for i, row in enumerate(rows) if row["converged"] == "T"] == [len(rows)-1]
    for row in rows[-2:]:
        assert row["orbital_ok"] == row["energy_ok"] == "T"
        assert float(row["max_scnsty"]) <= float(row["accy"])
        assert float(row["relative_energy_change"]) < .001 * float(row["accy"])
    for iteration in range(len(rows) + 1):
        assert (root / f"rwfn.out.iter{iteration:06d}").is_file()
    compare_wavefunctions(root / "rwfn.out", root / f"rwfn.out.iter{len(rows):06d}")
    print(f"{args.ranks} MPI rank(s): matched first sweep, restart converged in {len(rows)} cycles")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--tools", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path, required=True)
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--ranks", type=int, default=1)
    args = parser.parse_args()
    # The copied upstream memory kernels require Linux /proc/meminfo.
    # Do not silently substitute a different solver on macOS/Windows.
    if not Path('/proc/meminfo').is_file():
        print('SKIP: upstream memory MPI solver requires Linux /proc/meminfo')
        return 77
    for field in ("strict", "baseline", "tools", "mpiexec"):
        setattr(args, field, getattr(args, field).resolve())
    # Legacy MPI path helpers have short buffers. Use a short POSIX path.
    with tempfile.TemporaryDirectory(prefix="strict-mem-", dir="/tmp") as temporary:
        exercise(args, Path(temporary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
