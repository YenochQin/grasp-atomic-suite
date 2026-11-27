"""Real MCP producers and RMCDHF entry points; no SCF/physical claim."""

from __future__ import annotations

import argparse
import os
import struct
import subprocess
import tempfile
from pathlib import Path

from test_csf_descriptors import basic_blocks, write_fixture


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--angular", type=Path, required=True)
    parser.add_argument("--loader", type=Path, required=True)
    parser.add_argument("--mpiexec")
    args = parser.parse_args()
    angular, loader = args.angular.resolve(), args.loader.resolve()
    for ranks in [1, 2, 4] if args.mpiexec else [1]:
        with tempfile.TemporaryDirectory(prefix="rmcsf-") as temporary:
            work = Path(temporary)
            peel = ["2s", "2p-", "2p", "3s", "3p-", "3p"]
            parquet, header, text = write_fixture(work, ["1s"], peel, basic_blocks())
            parquet.rename(work / "rcsf_desc.parquet")
            header.rename(work / "rcsf_header.toml")
            text.unlink()
            environment = os.environ | {
                "MPI_TMP": str(work),
                "OPENBLAS_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1",
                "OMPI_MCA_rmaps_base_oversubscribe": "1",
            }
            prefix = [args.mpiexec, "-n", str(ranks)] if args.mpiexec else []
            result = subprocess.run(
                [*prefix, str(angular)],
                cwd=work,
                input="rcsf\nn\n" if args.mpiexec else "rcsf\ny\n",
                check=False,
                capture_output=True,
                text=True,
                env=environment,
                timeout=45,
            )
            assert result.returncode == 0 and "Execution complete" in result.stdout, (
                result.stdout + result.stderr
            )
            result = subprocess.run(
                [*prefix, str(loader)],
                cwd=work,
                check=False,
                capture_output=True,
                text=True,
                env=environment,
                timeout=30,
            )
            assert result.returncode == 0 and "PASS RMCDHF" in result.stdout, (
                result.stdout + result.stderr
            )
            expected = bytearray()
            for block in basic_blocks():
                for row in block:
                    values = {shell.label: shell.occupation for shell in row.shells}
                    expected.extend([2, *[values.get(label, 0) for label in peel]])
            for rank in range(ranks):
                assert (work / f"occupations.{rank:03d}.bin").read_bytes() == expected
            if ranks > 1:
                # Alter the electron count on just one rank.
                path = work / "001/mcp.mcpidx"
                original = path.read_bytes()
                data = bytearray(original)
                struct.pack_into("=q", data, 8 + 8*4 + 5*8, 99)
                path.write_bytes(data)
                result = subprocess.run(
                    [*prefix, str(loader)],
                    cwd=work,
                    check=False,
                    capture_output=True,
                    text=True,
                    env=environment,
                    timeout=30,
                )
                assert result.returncode != 0, "worker MCP mismatch accepted"
                path.write_bytes(original)
            header = work / "rcsf_header.toml"
            header.write_bytes(header.read_bytes() + b"\n")
            result = subprocess.run(
                [*prefix, str(loader)],
                cwd=work,
                check=False,
                capture_output=True,
                text=True,
                env=environment,
                timeout=30,
            )
            assert result.returncode != 0, "root header mismatch accepted"
    print("PASS RMCDHF native inputs, exact occupations and controlled errors")


if __name__ == "__main__":
    main()
