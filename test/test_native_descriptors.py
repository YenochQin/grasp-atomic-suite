"""Native Parquet reader integration, using the frozen text reader as oracle."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from test_csf_descriptors import Record, Shell, basic_blocks, write_fixture


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
        legacy_parquet = root / "csfs_descriptors.parquet"
        shutil.copyfile(parquet, legacy_parquet)
        rejected = subprocess.run(
            [binaries / "rangular", "--csf-descriptors", str(legacy_parquet), "--header", str(header)],
            cwd=root, input="y\n", capture_output=True, text=True, check=False, timeout=30,
        )
        assert rejected.returncode != 0 and "use a *_desc.parquet" in rejected.stderr, rejected
        shutil.copyfile(parquet, root / "rcsf_descriptors.parquet")
        shutil.copyfile(header, root / "rcsf_header.toml")
        rejected = subprocess.run(
            [binaries / "rangular"], cwd=root, input="y\n",
            capture_output=True, text=True, check=False, timeout=30,
        )
        assert rejected.returncode != 0 and "rcsf_desc.parquet" in rejected.stderr, rejected
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
        pair_parquet, _, pair_text = write_fixture(
            root / "pairs", [], ["1s", "2s", "3s", "4s"],
            [[
                Record((Shell("1s", 2), Shell("2s", 1, 1, 1), Shell("3s", 1, 1, 1, 0), Shell("4s", 1, 1, 1)), 1),
                Record((Shell("1s", 2), Shell("2s", 1, 1, 1), Shell("3s", 1, 1, 1, 2), Shell("4s", 1, 1, 1)), 1),
                Record((Shell("1s", 2), Shell("2s", 2), Shell("3s", 1, 1, 1)), 1),
            ]],
        )
        pair_text.unlink()
        result = subprocess.run(
            [str(lifecycle), str(pair_parquet), "pairs"],
            cwd=root, check=False, capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0 and "PASS descriptor lifecycle" in result.stdout, result.stdout + result.stderr
    print("Native descriptor integration passed")


if __name__ == "__main__":
    main()
