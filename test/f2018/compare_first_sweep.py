"""Compare fixed one-sweep synthetic runs with a separately built baseline.

Both binary directories must use the same compiler/math-library configuration.
This is a software-equivalence check, not an SCF-convergence calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "rmcdhf"))
from strict_mem_restart import inputs


ISODATA = (
    "Atomic number:\n10\nMass number (integer) :\n20\n"
    "Fermi distribution parameter a:\n0.52338755531043146\n"
    "Fermi distribution parameter c:\n3.0\nMass of nucleus (in amu):\n0\n"
    "Nuclear spin:\n0\nNuclear dipole moment:\n0\nNuclear quadrupole moment:\n0\n"
)
CSF = "Core subshells:\n\nPeel subshells:\n  1s\nCSF(s):\n  1s ( 2)\n\n         0+\n"


def run(binary: Path, stdin: str, work: Path) -> None:
    env = os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    (work / f"{binary.name}.stdin").write_text(stdin)
    result = subprocess.run([str(binary)], input=stdin, text=True, capture_output=True,
                            cwd=work, env=env, timeout=90)
    (work / f"{binary.name}.stdout.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{binary.name} exited {result.returncode}: {work}")
    if binary.name.startswith("rmcdhf"):
        assert re.findall(r"Iteration number\s+(\d+)", result.stdout) == ["1"], work


def exercise(bindir: Path, solver: str, work: Path) -> dict[str, str]:
    work.mkdir(parents=True)
    (work / "isodata").write_text(ISODATA)
    (work / "rcsf.inp").write_text(CSF)
    run(bindir / "rwfnestimate", "y\n2\n*\n", work)
    run(bindir / "rangular", "y\n", work)
    run(bindir / solver, inputs(1), work)
    files = [work / "rwfn.inp", work / "rwfn.out", work / "rmix.out", *sorted(work.glob("mcp.*"))]
    assert len(files) >= 6, "missing MCP records"
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-bindir", type=Path, required=True)
    parser.add_argument("--candidate-bindir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference = args.reference_bindir.resolve()
    candidate = args.candidate_bindir.resolve()
    results = {}
    for solver in ["rmcdhf", "rmcdhf_mem"]:
        old = exercise(reference, solver, output / solver / "reference")
        new = exercise(candidate, solver, output / solver / "candidate")
        results[solver] = {"reference": old, "candidate": new, "equal": old == new}
        (output / "comparison.json").write_text(json.dumps(results, indent=2) + "\n")
        assert old == new, f"one-sweep output changed for {solver}: {output}"
        print(f"{solver}: all {len(old)} wavefunction/mixing/MCP files are byte-identical")


if __name__ == "__main__":
    main()
