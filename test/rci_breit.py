"""Check the MPI Breit path with seven J=0+ CSFs from the suite's Ni I fixture."""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from native_fixture import descriptors
from mixing_stream import read_stream


def run_case(
    command: list[str], work: Path, inputs: Path, breit: bool
) -> float:
    shutil.copytree(inputs, work)
    # Keep the production frequency scale when the transverse term is enabled.
    transverse = "y\ny\n1e-6" if breit else "n\nn"
    result = subprocess.run(
        command,
        cwd=work,
        input=f"y\ncase\nn\n{transverse}\nn\nn\nn\nn\n1\n",
        capture_output=True,
        text=True,
        timeout=45,
        check=False,
        env=os.environ | {
            "MPI_TMP": str(work), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        },
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0 and "RCI_MPI: Execution complete." in output, output
    assert ("Breit integrals of type 1" in output) == breit, output
    # For ranks 1 through 7 the fixture admits 96, 97, 68, 34, 12, 3, 0
    # ordered orbital pairs. Type 2 combines two such pairs of equal rank.
    type2_counts = re.findall(r"Computing\s+(\d+)\s+Breit integrals of type 2", output)
    assert type2_counts == (["24558"] if breit else []), output
    block, = read_stream(work / "case.cimix", work / "case_desc.parquet")
    energy, = block.energies
    assert math.isfinite(energy), output
    return energy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rci", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path, required=True)
    parser.add_argument("--numproc-flag", required=True)
    args = parser.parse_args()
    command = [str(args.mpiexec.resolve()), args.numproc_flag, "2", str(args.rci.resolve())]
    data = Path(__file__).resolve().parents[1] / "data"
    with tempfile.TemporaryDirectory(prefix="rci-breit-") as temporary:
        root = Path(temporary)
        inputs = root / "inputs"
        inputs.mkdir()
        shutil.copyfile(data / "isodata", inputs / "isodata")
        shutil.copyfile(data / "test.wfn", inputs / "case.wfn")
        text = inputs / "case.c"
        with (data / "test.c").open() as source:
            text.write_text("".join(next(source) for _ in range(5 + 3 * 7)))
        descriptors(text)
        text.unlink()

        breit = run_case(command, root / "breit", inputs, breit=True)
        coulomb = run_case(command, root / "coulomb", inputs, breit=False)
    assert not math.isclose(breit, coulomb, rel_tol=1e-12), (breit, coulomb)
    print(
        f"Ni I J=0+ energies (hartree): Breit={breit:.12f}, "
        f"Coulomb={coulomb:.12f}, shift={breit - coulomb:+.12f}"
    )


if __name__ == "__main__":
    main()
