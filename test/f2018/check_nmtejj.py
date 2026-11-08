"""Compare the production lookup against a frozen pre-migration result."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main() -> None:
    reference = json.loads(Path(__file__).with_name("nmtejj_reference.json").read_text())
    result = subprocess.run([sys.argv[1]], capture_output=True, check=True, timeout=30)
    assert len(result.stdout.splitlines()) == reference["rows"]
    digest = hashlib.sha256(result.stdout).hexdigest()
    assert digest == reference["stdout_sha256"], f"NMTEJJ lookup changed: {digest}"
    print(f"NMTEJJ: {reference['rows']} cases match the legacy result exactly")


if __name__ == "__main__":
    main()
