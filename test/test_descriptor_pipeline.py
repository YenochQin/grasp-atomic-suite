"""Small biotransform/transition pipeline with generated one-electron fixtures."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from native_fixture import descriptors
from rdensity.regression import fixture


def execute(binary: Path, directory: Path, stdin: str, launcher: list[str]) -> None:
    result = subprocess.run(
        [*launcher, str(binary)],
        cwd=directory,
        input=stdin,
        text=True,
        check=False,
        capture_output=True,
        env=os.environ
        | {
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MPI_TMP": str(directory),
            "OMPI_MCA_rmaps_base_oversubscribe": "1",
        },
        timeout=60,
    )
    (directory / (binary.name + ".stdout")).write_text(result.stdout)
    (directory / (binary.name + ".stderr")).write_text(result.stderr)
    assert (
        result.returncode == 0 and "Execution complete" in result.stdout + result.stderr
    ), result.stdout + result.stderr


def pipeline(
    binaries: Path, root: Path, native: bool, launcher: list[str]
) -> dict[str, bytes]:
    root.mkdir(parents=True)
    template = root / "template"
    fixture(template, False, True)
    work = root / "biotransform"
    work.mkdir()
    shutil.copyfile(template / "isodata", work / "isodata")
    for stem in ("left", "right"):
        for extension in ("c", "w", "m"):
            shutil.copyfile(
                template / f"case.{extension}", work / f"{stem}.{extension}"
            )
        if native:
            descriptors(work / f"{stem}.c")
            (work / f"{stem}.c").unlink()
    suffix = "_mpi" if launcher else ""
    execute(
        binaries / ("rbiotransform" + suffix), work, "y\nn\nleft\nright\ny\n", launcher
    )
    outputs = {
        file.name: file.read_bytes()
        for pattern in ("*.bm", "*.bw")
        for file in work.glob(pattern)
    }
    assert len(outputs) == 4, outputs.keys()
    # RCI's existing solver requires at least as many CSFs in every block
    # as MPI ranks. Use a single two-CSF block for the two-rank regression.
    rci_work = root / "rci"
    fixture(rci_work, False, False)
    if native:
        descriptors(rci_work / "case.c")
        (rci_work / "case.c").unlink()
    execute(
        binaries / ("rci" + suffix),
        rci_work,
        "y\ncase\n" + ("n\n" if launcher and native else "") + "n\nn\nn\nn\nn\nn\n1\n",
        launcher,
    )
    outputs["rci.cm"] = (rci_work / "case.cm").read_bytes()
    for program, extensions in [
        ("gj", ["gj"]),
        ("rhfs", ["h", "hoffd"]),
        ("ris4", ["i"]),
    ]:
        answers = "y\ncase\nn\n" + ("n\nn\n" if program == "ris4" else "")
        execute(binaries / (program + suffix), rci_work, answers, launcher)
        for extension in extensions:
            path = rci_work / f"case.{extension}"
            assert path.stat().st_size > 0
            outputs[program + "." + extension] = path.read_bytes()
    for program in (
        ("rtransition_mpi",) if launcher else ("rtransition", "rtransition_phase")
    ):
        destination = root / program
        shutil.copytree(work, destination)
        execute(binaries / program, destination, "y\nn\nleft\nright\nE1\n", launcher)
        assert (destination / "left.right.t").stat().st_size > 0
        if native:
            assert not (destination / "SLASK").exists()
        outputs[program + ".t"] = (destination / "left.right.t").read_bytes()
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--mpiexec")
    parser.add_argument("--ranks", type=int, default=2)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="descriptor-pipeline-") as temporary:
        root = args.output.resolve() if args.output else Path(temporary)
        root.mkdir(exist_ok=True)
        launcher = [args.mpiexec, "-n", str(args.ranks)] if args.mpiexec else []
        candidate = pipeline(args.bin.resolve(), root / "native", True, launcher)
        if args.baseline:
            reference = pipeline(
                args.baseline.resolve(), root / "baseline", False, launcher
            )
            assert candidate == reference, "Pipeline output differs from baseline"
    print("PASS descriptor-only biotransform/transition pipeline")


if __name__ == "__main__":
    main()
