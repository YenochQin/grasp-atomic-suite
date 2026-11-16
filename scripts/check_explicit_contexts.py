"""Recheck built first-party Fortran contexts with implicit typing disabled.

Requires an already completed CMake build. Module outputs go to disposable
directories; the production module files and compiler standard are preserved.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shlex
import subprocess
import tempfile


def check(entry: dict[str, str]) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="explicit-context-") as temporary:
        arguments = shlex.split(entry["command"])
        command = []
        skip = False
        for argument in arguments:
            if skip:
                skip = False
                continue
            if argument in ("-o", "-J"):
                skip = True
            elif argument == "-c" or argument.startswith("-J"):
                continue
            else:
                command.append(argument)
        command += ["-fsyntax-only", "-fimplicit-none", "-J" + temporary]
        result = subprocess.run(command, cwd=entry["directory"], capture_output=True, text=True, timeout=60)
        return {"source": entry["file"], "object": entry["output"], "status": result.returncode,
                "diagnostics": result.stderr if result.returncode else ""}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    entries = json.loads((args.build / "compile_commands.json").read_text())
    entries = [entry for entry in entries if Path(entry["file"]).suffix.lower() == ".f90"
               and any(Path(entry["file"]).is_relative_to(root / name) for name in ("src", "test"))]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(check, entries))
    failures = [result for result in results if result["status"]]
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(f"{len(results)} contexts, {len({r['source'] for r in results})} sources, {len(failures)} failures")
    for result in failures:
        print(result["source"], result["diagnostics"])
    raise SystemExit(bool(failures))


if __name__ == "__main__":
    main()
