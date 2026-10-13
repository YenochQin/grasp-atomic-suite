"""Check production trace rows at a real ORTHY projection boundary."""
from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import tempfile
from pathlib import Path


def exercise(executable: Path, root: Path, option: str, rank: int) -> None:
    output = root / f"trace-{option}-{rank}"
    output.mkdir()
    result = subprocess.run(
        [str(executable), str(output), str(rank)],
        env=os.environ | {"GRASP_P0_TRACE": option},
        text=True, capture_output=True, check=True,
    )
    assert "checks passed" in result.stdout
    path = output / f"rmcdhf.p0.rank{rank}.csv"
    if option == "0" or option == "1" and rank != 0:
        assert not path.exists()
        return
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows and all(None not in row for row in rows)
    assert all(row["cycle"] == "2" and row["update"] == "1" for row in rows)
    assert all(row["rank"] == str(rank) and row["active_orbital"] == "2" for row in rows)
    stages = {row["stage"]: row for row in rows}
    assert float(stages["orth_before_projection"]["p0"]) == 2.0
    assert float(stages["orth_reference"]["p0"]) == 1.0
    projected = stages["orth_after_projection"]
    assert projected["peer"] == "1" and projected["inner_iteration"] == "3"
    assert float(projected["overlap"]) == 0.5
    assert float(projected["p0"]) == 1.5 and float(projected["p2"]) == 0.0
    normalized = stages["orth_normalized"]
    assert math.isclose(float(normalized["norm"]), 0.75, abs_tol=1e-14)
    assert math.isclose(float(normalized["scale"]), 1 / math.sqrt(0.75), abs_tol=1e-14)
    assert math.isclose(float(normalized["p0"]), math.sqrt(3), abs_tol=1e-14)
    assert math.isclose(float(normalized["p3"]), 1.0, abs_tol=1e-14)
    assert normalized["label"] == "2s"
    assert all(float(normalized[f"p{i}"]) == 0 for i in range(4, 8))
    assert math.isnan(float(stages["nonfinite_probe"]["p0"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    executable = parser.parse_args().executable.resolve()
    with tempfile.TemporaryDirectory(prefix="p0-trace-") as temporary:
        root = Path(temporary)
        for option, rank in (("0", 0), ("1", 1), ("1", 0), ("all", 1)):
            exercise(executable, root, option, rank)
        failure = subprocess.run(
            [str(executable), str(root / "missing" / "directory"), "0"],
            env=os.environ | {"GRASP_P0_TRACE": "1"}, text=True, capture_output=True, check=True,
        )
        assert "[P0_TRACE] cannot open" in failure.stderr
    print("P0 trace projection, rank selection, nonfinite values and disabled-output checks passed.")


if __name__ == "__main__":
    main()
