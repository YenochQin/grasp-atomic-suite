"""Native Parquet reader integration, using the frozen text reader as oracle."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from test_csf_descriptors import Record, basic_blocks, write_fixture


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, required=True)
    args = parser.parse_args()
    binaries = args.bin.resolve()
    with tempfile.TemporaryDirectory(prefix="csfn-") as temporary:
        root = Path(temporary)
        parquet, header, text = write_fixture(
            root, ["1s"], ["2s", "2p-", "2p", "3s", "3p-"], basic_blocks()
        )
        parquet = parquet.rename(root / "csfs_descriptors.parquet")
        result = subprocess.run(
            [str(binaries / "csf_descriptor_verify"), str(parquet), str(text)],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert "PASS" in result.stdout or "verified" in result.stdout.lower(), (
            result.stdout
        )
        lifecycle = binaries.parent / "test" / "test.descriptor_lifecycle"
        result = subprocess.run(
            [str(lifecycle), str(parquet)],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            env=os.environ | {"PATH": "/nonexistent"},
        )
        assert (
            result.returncode == 0 and "PASS descriptor lifecycle" in result.stdout
        ), result.stdout + result.stderr
        text.unlink()
        result = subprocess.run(
            [str(binaries / "rangular"), "--csf-descriptors", str(parquet)],
            cwd=root,
            input="y\n",
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            env=os.environ | {"PATH": "/nonexistent", "OPENBLAS_NUM_THREADS": "1"},
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert (root / "mcp.30").is_file(), result.stdout
        assert (root / "mcp.31").stat().st_size > 0
        # A valid legacy input must never rescue an invalid descriptor/header.
        write_fixture(
            root / "legacy", ["1s"], ["2s", "2p-", "2p", "3s", "3p-"], basic_blocks()
        )
        shutil.copyfile(root / "legacy" / "rcsf.inp", text)
        header.write_bytes(header.read_bytes() + b"\n")
        result = subprocess.run(
            [str(binaries / "rangular"), "--csf-descriptors", str(parquet)],
            cwd=root,
            input="y\n",
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode != 0, "header mismatch was accepted"
        core_parquet, _, core_text = write_fixture(
            root / "closed", ["1s"], [], [[Record((), 0)]]
        )
        core_text.unlink()
        result = subprocess.run(
            [str(lifecycle), str(core_parquet)],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            env=os.environ | {"PATH": "/nonexistent"},
        )
        assert (
            result.returncode == 0 and "PASS descriptor lifecycle" in result.stdout
        ), result.stdout + result.stderr
    print("Native descriptor integration passed")


if __name__ == "__main__":
    main()
