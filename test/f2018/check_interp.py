"""Production interpolation boundary/oracle checks, including a rejected empty table."""
from pathlib import Path
import subprocess
import sys


def main(binary: Path) -> None:
    for scenario in ("regular", "left", "short", "zero", "single", "single_node", "legacy", "empty"):
        result = subprocess.run([str(binary.resolve()), scenario], capture_output=True, text=True, timeout=10)
        if scenario == "empty":
            if result.returncode == 0 or "INTERP: empty table" not in result.stderr:
                raise RuntimeError(f"Wrong empty-table failure: {result.stdout}\n{result.stderr}")
        elif result.returncode or f"PASS interpolation {scenario}" not in result.stdout:
            raise RuntimeError(f"{scenario}: exit={result.returncode}\n{result.stdout}\n{result.stderr}")
        if scenario == "single" and "accuracy cannot be estimated" not in result.stdout:
            raise RuntimeError("One-point extrapolation must not claim that ACCY was verified")
        if scenario == "single_node" and "accuracy cannot be estimated" in result.stdout:
            raise RuntimeError("Exact table lookup must not be treated as extrapolation")
    print("PASS polynomial, first node, short table, zero, singleton, 71 regular points and rejected empty table")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
