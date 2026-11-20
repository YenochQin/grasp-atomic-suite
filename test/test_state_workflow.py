"""Exercise the state-named orbital, angular, SCF and CI programs in one directory."""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

from native_fixture import descriptors

CSF = """Core subshells:

Peel subshells:
  1s   2s
CSF(s):
  1s ( 2)

         0+
  2s ( 2)

         0+
"""


def workflow(binaries: Path, root: Path, launcher: list[str]) -> None:
    permanent = root / "serial io"
    permanent.mkdir()
    launch = root / "launch" if launcher else permanent
    scratch = root / "scratch"
    if launcher:
        launch.mkdir()
        scratch.mkdir()
        (launch / "disks").write_text(f"'{permanent}'\n'{scratch}'\n'{scratch}'\n")
    (permanent / "isodata").write_bytes(
        (Path(__file__).resolve().parents[1] / "data/isodata").read_bytes()
    )

    def run(program: str, answers: str, directory: Path | None = None) -> None:
        mpi = program.endswith("_mpi")
        result = subprocess.run(
            [*(launcher if mpi else []), binaries / program],
            cwd=directory or (launch if mpi else permanent),
            env=os.environ | {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
            input=answers,
            text=True,
            capture_output=True,
            timeout=45,
            check=False,
        )
        output = result.stdout + result.stderr
        assert result.returncode == 0 and "Execution complete." in output, output
        assert output.count("Name of state:") == 1, output

    suffix = "_mpi" if launcher else ""
    previous_files: dict[Path, bytes] = {}
    first_mcp: dict[str, bytes] = {}
    for name, solver in (("N", "rmcdhf"), ("N_3s", "rmcdhf_mem")):
        text = permanent / f"{name}.c"
        text.write_text(CSF if name == "N" else CSF.replace("2s", "3s"))
        descriptors(text)
        text.unlink()
        pair = [
            permanent / f"{name}{ending}"
            for ending in ("_header.toml", "_desc.parquet")
        ]
        if name == "N":
            older = permanent / f"{name}_descriptors.parquet"
            older.write_bytes(b"older descriptor must not shadow the current file")
            pair.append(older)
        previous_files.update({path: path.read_bytes() for path in pair})
        run("rwfnestimate", f"{name}\ny\n2\n*\n")
        run("rangular" + suffix, name + ("\nn\n" if launcher else "\ny\n"))
        run(solver + suffix, f"{name}\ny\n1\n1s\n\n1\n")
        run(
            "rci" + suffix,
            f"y\n{name}\n" + ("n\n" if launcher else "") + "n\n" * 6 + "1\n",
        )

        for path, expected in previous_files.items():
            assert path.read_bytes() == expected, path
        for extension in (
            ".initial.wfn",
            ".wfn",
            ".mix",
            ".rmcdhf.summary",
            ".cimix",
        ):
            path = permanent / (name + extension)
            assert path.is_file() and path.stat().st_size > 0, path
            previous_files[path] = path.read_bytes()
        for extension in (".rangular.log", ".rmcdhf.log"):
            assert (launch / (name + extension)).read_text().splitlines()[0] == name
        for rank in range(2 if launcher else 1):
            channel_dir = scratch / f"{rank:03d}" if launcher else permanent
            assert (channel_dir / "mcp.30").stat().st_size > 0
            assert not list(channel_dir.glob(f"{name}.mcp*"))
        if not launcher and name == "N":
            first_mcp = {p.name: p.read_bytes() for p in permanent.glob("mcp.*")}

    assert (permanent / "N.initial.wfn").read_bytes() != (
        permanent / "N_3s.initial.wfn"
    ).read_bytes()
    for generic in (
        "rcsf_header.toml",
        "rcsf_desc.parquet",
        "rcsf_descriptors.parquet",
        "rwfn.inp",
        "rwfn.out",
        "rmix.out",
        "rmcdhf.sum",
        "rmcdhf.log",
        "rangular.log",
        "mcp000.30",
    ):
        assert not list(root.rglob(generic)), generic

    if not launcher:
        (permanent / "N.initial_header.toml").write_bytes(
            (permanent / "N_header.toml").read_bytes()
        )
        path = permanent / "N.initial.wfn"
        original = path.read_bytes()
        for program, state in (("rwfnestimate", "N"), ("rmcdhf", "N.initial")):
            result = subprocess.run(
                [binaries / program],
                cwd=permanent,
                input=f"{state}\ny\n",
                text=True,
                capture_output=True,
                timeout=15,
                check=False,
            )
            assert (
                result.returncode != 0 and "Orbital path also belongs" in result.stderr
            ), result.stderr
            assert path.read_bytes() == original
        # The same calculation with a legacy answer file must produce identical binary records.
        legacy = root / "legacy"
        legacy.mkdir()
        for source, target in (
            ("N_header.toml", "rcsf_header.toml"),
            ("N_desc.parquet", "rcsf_descriptors.parquet"),
            ("isodata", "isodata"),
        ):
            (legacy / target).write_bytes((permanent / source).read_bytes())
        for program, answers in (
            ("rwfnestimate", "y\n2\n*\n"),
            ("rangular", "y\n"),
            ("rmcdhf", "y\n1\n1s\n\n1\n"),
        ):
            run(program, answers, legacy)
        for source, target in (
            ("rwfn.inp", "N.initial.wfn"),
            ("rwfn.out", "N.wfn"),
            ("rmix.out", "N.mix"),
        ):
            assert (legacy / source).read_bytes() == (permanent / target).read_bytes()
        assert {p.name: p.read_bytes() for p in legacy.glob("mcp.*")} == first_mcp


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("binaries", type=Path)
    parser.add_argument("--mpiexec")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="states-") as directory:
        workflow(
            args.binaries.resolve(),
            Path(directory),
            [args.mpiexec, "-n", "2"] if args.mpiexec else [],
        )
