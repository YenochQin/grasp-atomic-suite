"""Native Parquet reader integration, using the frozen text reader as oracle."""

from __future__ import annotations

import argparse
from itertools import combinations, islice
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
        alternate, _, _ = write_fixture(
            root / "alternate", ["1s"], [], [[Record((), 0)]]
        )
        result = subprocess.run(
            [str(lifecycle), str(parquet), "dual", str(alternate)],
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
            input="case\ny\n",
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            env=os.environ | {"PATH": "/nonexistent", "OPENBLAS_NUM_THREADS": "1"},
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert (root / "case.mcp.30").is_file(), result.stdout
        assert (root / "case.mcp.31").stat().st_size > 0
        legacy_parquet = root / "csfs_descriptors.parquet"
        shutil.copyfile(parquet, legacy_parquet)
        rejected = subprocess.run(
            [binaries / "rangular", "--csf-descriptors", str(legacy_parquet), "--header", str(header)],
            cwd=root, input="case\ny\n", capture_output=True, text=True, check=False, timeout=30,
        )
        assert rejected.returncode != 0 and "use a *_desc.parquet" in rejected.stderr, rejected
        shutil.copyfile(parquet, root / "case_descriptors.parquet")
        shutil.copyfile(header, root / "case_header.toml")
        rejected = subprocess.run(
            [binaries / "rangular"], cwd=root, input="case\ny\n",
            capture_output=True, text=True, check=False, timeout=30,
        )
        assert rejected.returncode != 0 and "case_desc.parquet" in rejected.stderr, rejected
        # A valid legacy input must never rescue an invalid descriptor/header.
        write_fixture(
            root / "legacy", ["1s"], ["2s", "2p-", "2p", "3s", "3p-"], basic_blocks()
        )
        shutil.copyfile(root / "legacy" / "rcsf.inp", text)
        header.write_bytes(header.read_bytes() + b"\n")
        result = subprocess.run(
            [str(binaries / "rangular"), "--csf-descriptors", str(parquet)],
            cwd=root,
            input="case\ny\n",
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
        orbital_numbers = list(range(1, 17))
        batch = [
            Record(tuple(Shell(f"{n}s", 2) if n in paired else Shell(f"{n}s", 1, 1, 1)
                         for n in sorted(paired + single)), 0)
            for paired, single in islice(
                ((paired, single) for paired in combinations(orbital_numbers, 2)
                 for single in combinations([n for n in orbital_numbers if n not in paired], 2)),
                8193,
            )
        ]
        batch_parquet, _, _ = write_fixture(root / "batches", [], [f"{n}s" for n in orbital_numbers], [batch])
        result = subprocess.run(
            [str(lifecycle), str(batch_parquet), "batches"],
            cwd=root, check=False, capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0 and "PASS descriptor batch boundary" in result.stdout, result.stdout + result.stderr
    print("Native descriptor integration passed")


if __name__ == "__main__":
    main()
