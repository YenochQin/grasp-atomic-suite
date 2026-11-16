"""Use one existing CSF to distinguish normal EOF from a missing required record."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("executable", type=Path)
parser.add_argument("--mpiexec")
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
source = (root/"data/test.c").read_text().splitlines(keepends=True)
command = [str(args.executable.resolve())]
if args.mpiexec:
    command = [args.mpiexec, "-n", "2", *command]
with tempfile.TemporaryDirectory(prefix="csf-record-io-") as temporary:
    for case, length in (("valid", 8), ("missing-subshell", 6), ("missing-coupling", 7)):
        work = Path(temporary)/case
        work.mkdir()
        (work/"valid.c").write_text("".join(source[:8]))
        (work/"truncated.c").write_text("".join(source[:length]))
        result = subprocess.run([*command,case], cwd=work, capture_output=True, text=True, timeout=15)
        output = result.stdout+result.stderr
        if case == "valid":
            assert result.returncode == 0 and "PASS CSF record input" in output, output
        else:
            assert result.returncode != 0 and "Input I/O failure: CSF" in output and "unit=21" in output, output
        print(case + (": accepted" if case == "valid" else ": rejected before using missing data"))
