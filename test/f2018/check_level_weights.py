"""Exercise the production weight-input branch, including a retry."""
from pathlib import Path
import subprocess
import sys


def main(binary: Path) -> None:
    cases = {"equal": "1\n", "standard": "5\n", "custom": "9\n2 3 5\n",
             "retry": "9\n-1 3 5\n2 3 5\n"}
    for scenario, answers in cases.items():
        result = subprocess.run([str(binary.resolve()), scenario], input=answers, text=True,
                                capture_output=True, timeout=10)
        if result.returncode or f"PASS level weights {scenario}" not in result.stdout:
            raise RuntimeError(f"{scenario}: exit={result.returncode}\n{result.stdout}\n{result.stderr}")
    print("PASS equal, standard, custom and retry weights; exact extent and normalized values")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
