"""Check the valid Breit coefficients and diagnostics for invalid inputs."""
import subprocess
import sys

command = sys.argv[1:]
if len(command) == 1:
    subprocess.run([*command, "valid"], timeout=15, check=True)
for scenario, diagnostic in [
    ("orbital", "invalid orbital combination"),
    ("rank-zero", "Invalid tensor_rank = 0"),
]:
    result = subprocess.run([*command, scenario], timeout=15, capture_output=True, text=True)
    assert result.returncode != 0 and diagnostic in result.stderr, result.stdout + result.stderr
    assert "Accepted invalid Breit input" not in result.stderr, result.stderr
