#!/usr/bin/env python3
"""Run the native array comparison, rejecting premature legacy STOP exits."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SUCCESS_MARKER = "CSF_DESCRIPTOR_VERIFICATION_SUCCESS"


def verify_descriptors(
    verifier: Path, binary: Path, text: Path
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [str(verifier.resolve()), str(binary.resolve()), str(text.resolve())],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode or SUCCESS_MARKER not in result.stdout.splitlines():
        raise ValueError(
            f"CSF array verification failed (exit {result.returncode}):\n{result.stdout}{result.stderr}"
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("text", type=Path)
    parser.add_argument("--verifier", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        result = verify_descriptors(
            arguments.verifier, arguments.binary, arguments.text
        )
    except (OSError, ValueError) as error:
        print(error, file=sys.stderr)
        return 1
    print(result.stdout, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
