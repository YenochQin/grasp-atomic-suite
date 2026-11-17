"""Exercise orbital insertion and independent left/right block boundaries."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

from test_csf_descriptors import Record, Shell, write_fixture


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("--mpiexec")
    args = parser.parse_args()
    for scenario in ("single", "multi", "invalid"):
        with tempfile.TemporaryDirectory(prefix="csf-merge-") as temporary:
            work = Path(temporary)
            right_peel = ["2s", "4s", "2p"]
            right_blocks = [
                [Record((Shell("4s", 2),), 0)],
                [Record((Shell("2s", 1, 1, 1), Shell("2p", 1, 3, 1)), 2, -1)],
            ]
            if scenario == "multi":
                right_peel = ["2s", "4s", "5s", "2p"]
                right_blocks[1].insert(
                    0, Record((Shell("5s", 1, 1, 1), Shell("2p", 1, 3, 1)), 2, -1)
                )
            elif scenario == "invalid":
                right_peel = ["2p", "2s", "4s"]
                right_blocks[1] = [
                    Record((Shell("2p", 1, 3, 1), Shell("2s", 1, 1, 1)), 2, -1)
                ]
            for stem, peel, blocks in [
                (
                    "left",
                    ["2s", "3s", "2p"],
                    [
                        [Record((Shell("2s", 2),), 0)],
                        [Record((Shell("3s", 1, 1, 1), Shell("2p", 1, 3, 1)), 2, -1)],
                    ],
                ),
                ("right", right_peel, right_blocks),
            ]:
                parquet, header, text = write_fixture(work / stem, ["1s"], peel, blocks)
                parquet.rename(work / f"{stem}_descriptors.parquet")
                header.rename(work / f"{stem}_header.toml")
                text.unlink()
            command = [str(args.executable.resolve()), scenario]
            if args.mpiexec:
                command = [args.mpiexec, "-n", "2", *command]
            result = subprocess.run(
                command,
                cwd=work,
                check=False,
                capture_output=True,
                text=True,
                timeout=30,
            )
            if scenario == "invalid":
                assert (
                    result.returncode != 0
                    and "inconsistent orbital order" in result.stdout + result.stderr
                ), result.stdout + result.stderr
            else:
                assert result.returncode == 0 and "PASS unequal" in result.stdout, (
                    result.stdout + result.stderr
                )
            assert not (work / "SLASK").exists()
    print("PASS dual descriptor mapping")


if __name__ == "__main__":
    main()
