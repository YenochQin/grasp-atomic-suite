"""Verify START tracing against the unchanged library routine and CSV values."""
from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import tempfile
from pathlib import Path


def exercise(executable: Path, root: Path, mode: str, rank: int) -> None:
    directory = root / f"{mode}-{rank}"
    directory.mkdir()
    result = subprocess.run(
        [str(executable), str(directory), str(rank)],
        env=os.environ | {"GRASP_START_TRACE": mode},
        capture_output=True, text=True, check=True,
    )
    assert "START parity checks passed" in result.stdout
    path = directory / f"rmcdhf.start.rank{rank}.csv"
    if mode == "0" or mode == "1" and rank != 0:
        assert not list(directory.iterdir())
        return
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    with (directory / f"rmcdhf.start_terms.rank{rank}.csv").open() as stream:
        terms = list(csv.DictReader(stream))
    assert all(None not in row for row in rows + terms)
    assert {r["itype"] for r in rows} == {"1", "2", "3"}
    assert {r["nparm"] for r in rows} == {"0", "2"}
    assert {r["label"] for r in rows} == {"3d", "3d-"}
    assert {r["rank"] for r in rows} == {str(rank)}
    assert {r["stage"] for r in rows} == {"initial", "iteration", "final"}
    final = [r for r in rows if r["stage"] == "final"]
    assert any(r["fallback"] == "T" and r["start_iteration"] == "36" for r in final)
    assert any(r["fallback"] == "F" and r["meets_accy"] == "T" for r in final)
    for row in rows:
        if row["point"] == "6" and row["stage"] != "initial":
            expected = float(row["delta_scaled"]) * float(row["r"]) ** float(row["gamma"])
            assert math.isclose(float(row["delta_physical"]), expected, rel_tol=1e-13, abs_tol=1e-300)
    assert {r["kind"] for r in terms} == {"cnc", "coupling"}
    coupling = [r for r in terms if r["kind"] == "coupling"]
    assert all(r["itype"] == "2" and r["peer"] == "1" for r in coupling)
    assert all(float(r["value"]) == 0.02 and float(r["peer_p0"]) == 221.0 for r in coupling)
    assert {int(r["peer"]) for r in terms if r["kind"] == "cnc"} == set(range(2, 7))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    executable = parser.parse_args().executable.resolve()
    with tempfile.TemporaryDirectory(prefix="start-trace-") as temporary:
        root = Path(temporary)
        for mode, rank in [("0", 0), ("1", 1), ("1", 0), ("all", 1)]:
            exercise(executable, root, mode, rank)
        result = subprocess.run(
            [str(executable), str(root / "missing"), "0"],
            env=os.environ | {"GRASP_START_TRACE": "1"},
            capture_output=True, text=True, check=True,
        )
        assert "[START_TRACE] cannot open" in result.stderr
        assert "START parity checks passed" in result.stdout
    print("START exact parity, rank controls, coefficients, iterations, fallback and I/O failure passed.")


if __name__ == "__main__":
    main()
