"""A single rank must reject an oversized count before entering a payload collective."""
import subprocess
import sys

result = subprocess.run(sys.argv[1:], capture_output=True, text=True, timeout=20, check=False)
expected = "MPI count probe on rank 1: element count outside default MPI range"
if result.returncode <= 0 or expected not in result.stderr:
    raise RuntimeError(f"Wrong failure (exit {result.returncode}):\n{result.stdout}\n{result.stderr}")
print(expected)
