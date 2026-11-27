"""Regress neutral Ni/Zn TF estimates, including first h and i orbitals.

Uses the default finite-nucleus grid and all orbitals through n=10, l=6.
rwfnestimate reads orbital metadata from the header TOML.
Requires Python 3.8+ and no third-party packages.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import json
import tempfile


def field(n: int, symbol: str) -> str:
    return f"{n:3d}{symbol:<2}"


def header(z: int) -> tuple[str, set[tuple[int, int]]]:
    core = [field(1, "s"), field(2, "s"), field(2, "p-"), field(2, "p"),
            field(3, "s"), field(3, "p-"), field(3, "p")]
    peel = []
    expected = set()
    for n in range(1, 11):
        for l, symbol in enumerate("spdfghi"[:min(n, 7)]):
            for kappa in ([l, -l - 1] if l else [-1]):
                label = field(n, symbol + ("-" if kappa > 0 else ""))
                expected.add((n, kappa))
                if label not in core:
                    peel.append(label)
    occupation = [("d-", 4), ("d", z - 24)]
    first = "".join(f"{field(3, symbol)}({q:2d})" for symbol, q in occupation)
    first += f"{field(4, 's')}( 2)"
    return "\n".join(["Core subshells:", "".join(core), "Peel subshells:",
                      "".join(peel), "CSF(s):", first]) + "\n", expected


def run(binary: Path, directory: Path, stdin: str) -> str:
    (directory / f"{binary.name}.stdin").write_text(stdin)
    result = subprocess.run([str(binary)], cwd=directory, input=stdin,
                            text=True, capture_output=True, timeout=30, check=False)
    log = result.stdout + result.stderr
    (directory / f"{binary.name}.log").write_text(log)
    if result.returncode:
        raise AssertionError(f"{binary.name} exited {result.returncode}\n{log}")
    return log


# Share the independent binary reader with RMCDHF regression utilities.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "rmcdhf"))
from compare_rwfn import read_radial_stream


def norm(p: tuple[float, ...], q: tuple[float, ...], r: tuple[float, ...],
         z: int) -> float:
    # Independent Simpson integration in x=log(1+r/RNT), with dr/dx=r+RNT.
    # Extend the negligible truncated tail with zeros for an even interval count.
    y = [(a * a + b * b) * (radius + 2e-6 / z)
         for a, b, radius in zip(p, q, r)]
    y.extend([0.0] * (3 if len(y) % 2 == 0 else 2))
    return 0.05 / 3 * (y[0] + y[-1] + 4 * sum(y[1:-1:2]) + 2 * sum(y[2:-1:2]))


def validate(path: Path, expected: set[tuple[int, int]], z: int) -> list[dict]:
    orbitals = read_radial_stream(path)
    summary = []
    for (n, kappa), orbital in orbitals.items():
        p, q, r, energy = orbital.p, orbital.q, orbital.r, orbital.energy
        assert all(math.isfinite(value) for value in (energy, *p, *q, *r)), (n, kappa)
        assert energy > 0 and any(p), (n, kappa)
        assert r[0] == 0 and all(a < b for a, b in zip(r, r[1:])), (n, kappa)
        integral = norm(p, q, r, z)
        # Allow for independent quadrature and the program's H**6 accuracy.
        assert abs(integral - 1) < 100 * 0.05**6, (n, kappa, integral)
        if (n, kappa) in {(6, 5), (6, -6), (7, 6), (7, -7)}:
            # These first h/i large components must be nodeless. Ignore
            # amplitudes below SOLVH's own final truncation accuracy.
            significant = [value for value in p if abs(value) > 0.1 * 0.05**6]
            assert significant and all(value * significant[0] > 0 for value in significant)
        summary.append({"n": n, "kappa": kappa, "energy": energy, "norm": integral})
    seen = set(orbitals)
    assert seen == expected, (seen - expected, expected - seen)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bindir", type=Path, required=True)
    parser.add_argument("--element", choices=["Ni", "Zn"], required=True)
    parser.add_argument("--workdir", type=Path)
    args = parser.parse_args()
    bindir = args.bindir.resolve()
    directory = args.workdir or Path(tempfile.mkdtemp(prefix=f"rwfnestimate-{args.element}-"))
    directory.mkdir(parents=True, exist_ok=True)
    print(f"Inputs, logs and orbitals: {directory}", flush=True)
    z = {"Ni": 28, "Zn": 30}[args.element]
    content, expected = header(z)
    (directory / "case_header.toml").write_text("[header_info]\nheader_lines = " + json.dumps(content.splitlines()[:5]) + "\nelectron_count = " + str(z) + "\n")
    # Static, spin-zero finite nuclei generated with upstream rnucleus.
    shutil.copyfile(Path(__file__).parent / "fixtures" / f"{args.element}.isodata",
                    directory / "isodata")
    log = run(bindir / "rwfnestimate", directory, "case\ny\n2\n*\n")
    assert "RWFNESTIMATE: Execution complete." in log
    for error in (
        "Unable to compute",
        "calculate_radial_dirac_coefficients: Join",
        "integrate_radial_tail: Grid",
        "could not be decoded",
    ):
        assert error not in log, log
    summary = validate(directory / "case.initial.wfn", expected, z)
    (directory / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Validated {len(summary)} Thomas-Fermi orbitals for neutral {args.element}.")


if __name__ == "__main__":
    main()
