"""Run real platform wrappers in a disposable directory."""
from __future__ import annotations

import subprocess
from pathlib import Path
import sys
import tempfile


with tempfile.TemporaryDirectory(prefix="suite-paths-") as temporary:
    subprocess.run([sys.argv[1]], cwd=temporary, check=True, timeout=30)
    exited = subprocess.run([sys.argv[1], "exit"], cwd=temporary, timeout=30)
    assert exited.returncode == 7, exited.returncode
    assert (Path(temporary) / "exit-flush").read_text() == "buffered output\n"
    exited = subprocess.run([sys.argv[1], "zero"], cwd=temporary, timeout=30)
    assert exited.returncode == 0, exited.returncode
