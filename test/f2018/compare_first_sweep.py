"""Compare fixed one-sweep runs, including complete records and diagnostics.

The supplied builds must have matching compiler, math library and options.
Binary hashes record identity; they do not independently certify provenance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "rmcdhf"))
from strict_mem_restart import inputs
from compare_rwfn import detect_endian, read_record, read_rwfn

ISODATA = (
    "Atomic number:\n10\nMass number (integer) :\n20\n"
    "Fermi distribution parameter a:\n0.52338755531043146\n"
    "Fermi distribution parameter c:\n3.0\nMass of nucleus (in amu):\n0\n"
    "Nuclear spin:\n0\nNuclear dipole moment:\n0\nNuclear quadrupole moment:\n0\n"
)
CSF = "Core subshells:\n\nPeel subshells:\n  1s\nCSF(s):\n  1s ( 2)\n\n         0+\n"
BINARY_FILES = {"rwfn.inp", "rwfn.out", "rmix.out", "mcp.30", "mcp.31", "mcp.32", "mcp.33"}
NATIVE_TEXT_FILES = {"rangular.log", "rmcdhf.log", "rmcdhf.sum"}
MCP_RECORD_SIZES = {
    "mcp.30": [12, 4, 8, 15, 12, 12, 11, 4, 8],
    "mcp.31": [15, 12, 12, 15, 8],
    "mcp.32": [15, 12, 12, 15, 8],
    "mcp.33": [15, 12, 12, 15, 8],
}
THREADS = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_DYNAMIC": "FALSE"}
# Only these anchored, complete metadata records may differ. Scientific lines,
# warnings, stopping messages, iteration counts and runtime diagnostics stay.
METADATA = (
    r"\s*Date \(Yr/Mon/Day\): \d{4}/\d{2}/\d{2}\s*",
    r"\s*Time \(Hr/Min/Sec\): \d{2}/\d{2}/\d{2}\.\d+\s*",
    r"\s*Zone: [+-]\d{4}\s*",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def normalize(text: str) -> str:
    result: list[str] = []
    wall_time = False
    for line in text.splitlines(keepends=True):
        if wall_time and re.fullmatch(r"\s*\d+ seconds\s*", line):
            result.append("<wall seconds>\n")
        elif any(re.fullmatch(pattern, line) for pattern in METADATA):
            result.append("<date/time metadata>\n")
        else:
            result.append(line)
        wall_time = line.strip() == "Wall time:"
    return "".join(result)


def check_diagnostics(text: str) -> None:
    if re.search(r"(?<![A-Za-z])(?:NaN|[+-]?Inf(?:inity)?)(?![A-Za-z])", text, re.I):
        raise ValueError("Nonfinite diagnostic")
    if re.search(r"Fortran runtime error|ERROR STOP|IEEE_(?:INVALID|ZERO|OVERFLOW)_FLAG", text, re.I):
        raise ValueError("Runtime failure diagnostic")


def run(binary: Path, stdin: str, work: Path) -> dict[str, object]:
    (work / f"{binary.name}.stdin").write_text(stdin)
    try:
        result = subprocess.run([str(binary)], input=stdin, text=True, capture_output=True,
                                cwd=work, env=os.environ | THREADS, timeout=90)
    except subprocess.TimeoutExpired as error:
        for name, data in (("stdout", error.stdout), ("stderr", error.stderr)):
            (work / f"{binary.name}.{name}.log").write_bytes(data.encode() if isinstance(data, str) else data or b"")
        raise RuntimeError(f"{binary.name} timed out: {work}") from error
    for name, data in (("stdout", result.stdout), ("stderr", result.stderr)):
        (work / f"{binary.name}.{name}.log").write_text(data)
    for data in (result.stdout, result.stderr):
        check_diagnostics(data)
    if result.returncode:
        raise RuntimeError(f"{binary.name} exited {result.returncode}: {work}")
    combined = result.stdout + result.stderr
    program = "RMCDHF" if binary.name.startswith("rmcdhf") else binary.name.upper()
    if f"{program}: Execution complete." not in combined:
        raise ValueError(f"{binary.name}: missing completion marker")
    if binary.name.startswith("rmcdhf"):
        if re.findall(r"Iteration number\s+(\d+)", combined) != ["1"]:
            raise ValueError("Expected exactly one SCF iteration")
        if combined.count("Maximum iterations in SCF Exceeded.") != 1:
            raise ValueError("Expected the fixed one-iteration limit, not a convergence claim")
    return {"returncode": result.returncode, "stdout": normalize(result.stdout), "stderr": normalize(result.stderr)}


def check_binary_files(work: Path) -> None:
    actual = {p.name for p in work.glob("mcp.*")}
    if actual != {"mcp.30", "mcp.31", "mcp.32", "mcp.33"}:
        raise ValueError(f"Incomplete/unexpected MCP file set: {sorted(actual)}")
    for name in BINARY_FILES:
        path = work / name
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing/empty required output: {name}")
        endian = detect_endian(path) if name.startswith(("rwfn", "rmix")) else "<"
        with path.open("rb") as stream:
            records = []
            while (record := read_record(stream, endian)) is not None:
                records.append(record)
        if name in MCP_RECORD_SIZES and [len(record) for record in records] != MCP_RECORD_SIZES[name]:
            raise ValueError(f"Incomplete/unexpected MCP records for the frozen single-CSF fixture: {name}")
        if name.startswith("rwfn"):
            orbitals = read_rwfn(path)
            if set(orbitals) != {(1, -1)}:
                raise ValueError(f"Unexpected orbital identity: {name}")
            orbital = orbitals[1, -1]
            if not orbital.r or not all(math.isfinite(v) for v in (orbital.energy, *orbital.r, *orbital.p, *orbital.q)):
                raise ValueError(f"Nonfinite/empty orbital: {name}")
        if name == "rmix.out":
            if len(records) != 6 or records[0] != b"G92MIX":
                raise ValueError("Incomplete mixing records")
            if struct.unpack(f"{endian}6i", records[1]) != (2, 1, 1, 1, 1, 1):
                raise ValueError("Wrong mixing identity/counts")
            if struct.unpack(f"{endian}5i", records[2]) != (1, 1, 1, 1, 1):
                raise ValueError("Wrong mixing block")
            if struct.unpack(f"{endian}i", records[3]) != (1,):
                raise ValueError("Wrong root identity")
            values = struct.unpack(f"{endian}2d", records[4]) + struct.unpack(f"{endian}d", records[5])
            if not all(map(math.isfinite, values)):
                raise ValueError("Nonfinite mixing coefficients/energy")


def exercise(bindir: Path, solver: str, work: Path) -> dict:
    work.mkdir(parents=True)
    (work / "isodata").write_text(ISODATA)
    (work / "rcsf.inp").write_text(CSF)
    runs = {}
    for binary, stdin in (("rwfnestimate", "y\n2\n*\n"), ("rangular", "y\n"), (solver, inputs(1))):
        runs[binary] = run(bindir / binary, stdin, work)
    check_binary_files(work)
    texts = {}
    for name in sorted(NATIVE_TEXT_FILES):
        path = work / name
        if not path.is_file():
            raise ValueError(f"Missing scientific text output: {name}")
        text = path.read_text()
        check_diagnostics(text)
        texts[name] = normalize(text)
    expected = BINARY_FILES | NATIVE_TEXT_FILES | {"isodata", "rcsf.inp"}
    for binary in runs:
        expected |= {f"{binary}.stdin", f"{binary}.stdout.log", f"{binary}.stderr.log"}
    actual = {p.name for p in work.iterdir()}
    if actual != expected:
        raise ValueError(f"Unexpected output inventory: missing={expected-actual}, extra={actual-expected}")
    return {"binaries": {name: sha256((work / name).read_bytes()) for name in sorted(BINARY_FILES)},
            "runs": runs, "scientific_text": texts}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-bindir", type=Path, required=True)
    parser.add_argument("--candidate-bindir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--solvers", nargs="+", choices=("rmcdhf", "rmcdhf_mem"), default=["rmcdhf", "rmcdhf_mem"],
                        help="Explicit partial scope is recorded; the default requires both solvers")
    args = parser.parse_args()
    output = args.output.resolve()
    reference, candidate = args.reference_bindir.resolve(), args.candidate_bindir.resolve()
    if reference == candidate:
        parser.error("Reference and candidate must be independent build directories")
    output.mkdir(parents=True, exist_ok=False)
    report: dict = {"status": "running", "solvers": args.solvers, "full_serial_scope": set(args.solvers)=={"rmcdhf","rmcdhf_mem"},
                    "threads": THREADS, "inputs": {"isodata": sha256(ISODATA.encode()), "rcsf.inp": sha256(CSF.encode()),
                    "scf_stdin": sha256(inputs(1).encode())}, "results": {}}
    try:
        report["executables"] = {
            side: {name: sha256((bindir/name).read_bytes()) for name in ["rwfnestimate", "rangular", *args.solvers]}
            for side, bindir in (("reference", reference), ("candidate", candidate))}
        for solver in args.solvers:
            old = exercise(reference, solver, output / solver / "reference")
            new = exercise(candidate, solver, output / solver / "candidate")
            report["results"][solver] = {"reference": old, "candidate": new, "equal": old == new}
            if old != new:
                raise ValueError(f"Scientific output/diagnostic difference for {solver}")
            print(f"{solver}: seven binary outputs and all scientific text/diagnostics match")
        report["status"] = "passed"
        return 0
    except (OSError, RuntimeError, ValueError, struct.error) as error:
        report.update(status="failed", error=str(error))
        print(str(error), file=sys.stderr)
        return 1
    finally:
        (output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    sys.exit(main())
