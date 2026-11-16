"""Require controlled failures on bad MCP input and exact valid/empty payloads."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("--mpiexec")
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--reader", action="store_true")
    args = parser.parse_args()
    command = [str(args.executable.resolve())]
    if args.mpiexec:
        command = [args.mpiexec, "-n", str(args.ranks), *command]
    cases = (("valid", "empty", "missing", "global-header", "block-sizes", "block-labels",
              "channel-header", "basis", "flags", "magic", "rank", "sparse-block-id", "block-id",
              "label", "payload", "sentinel", "short-global-header", "short-label", "short-payload")
             if args.reader else ("valid", "empty", "label", "sentinel", "block-header", "payload",
                                  "load-label", "load-sentinel", "load-block-header", "load-payload",
                                  "short-label", "short-payload", "load-short-label", "load-short-payload"))
    with tempfile.TemporaryDirectory(prefix="mcp-io-") as temporary:
        for case in cases:
            folder = Path(temporary) / case
            folder.mkdir()
            result = subprocess.run(
                [*command, case], cwd=folder, capture_output=True, text=True, timeout=15,
                env=os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
            )
            output = result.stdout + result.stderr
            if case in {"valid", "empty"}:
                marker = "PASS MCP reader input" if args.reader else "PASS MCP memory input"
                assert result.returncode == 0 and marker in output, (case, output)
            else:
                assert result.returncode != 0 and "Input I/O failure: MCP" in output, (case, output)
                if args.reader:
                    unit = 30 if case in {"missing", "global-header", "short-global-header", "block-sizes", "block-labels", "sparse-block-id"} else 31
                    assert f"unit={unit}" in output and any(name in output for name in ("setmcp", "setcof")), (case, output)
                else:
                    assert "setmcp_mem" in output and "unit=32" in output, (case, output)
                assert "Corrupt MCP was accepted" not in output, (case, output)
                if "short" in case:
                    status = re.search(r"iostat=(-?\d+)", output)
                    assert status and int(status[1]) > 0, (case, output)
            print(f"{case}: {'accepted' if case in {'valid', 'empty'} else 'rejected with context'}")


if __name__ == "__main__":
    main()
