"""Exercise the state-named orbital, angular, SCF and CI programs in one directory."""

import argparse
import os
import subprocess
import struct
import tempfile
from pathlib import Path

from native_fixture import descriptors
from rdensity.regression import records, write_records
from test_csf_descriptors import Record, Shell, text_record

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

    def run(
        program: str, answers: str, directory: Path | None = None, *, failure: str | None = None
    ) -> None:
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
        if failure:
            assert result.returncode != 0 and failure in output, output
            assert "Execution complete." not in output, output
        else:
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
            assert (channel_dir / ("mcp.mcpidx" if launcher else "mcp.30")).stat().st_size > 0
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

    saved_wave = permanent / "N_3s.wfn"
    archived_wave = saved_wave.with_suffix(".w")
    saved_wave.rename(archived_wave)
    result = subprocess.run(
        [*launcher, binaries / ("rci" + suffix)],
        cwd=launch,
        input="y\nN_3s\n" + ("n\n" if launcher else "") + "n\n" * 6 + "1\n",
        text=True,
        capture_output=True,
        timeout=45,
        check=False,
    )
    assert result.returncode != 0 and "N_3s.wfn" in result.stderr, result.stderr
    assert not saved_wave.exists()
    assert archived_wave.read_bytes() == previous_files[saved_wave]

    if not launcher:
        result = subprocess.run(
            [binaries / "rwfnestimate"],
            cwd=permanent,
            input="N_3s\ny\n1\nN_3s.w\n*\n",
            text=True,
            capture_output=True,
            timeout=15,
            check=False,
        )
        assert result.returncode != 0 and "Use a .wfn" in result.stderr, result.stderr
        assert (permanent / "N_3s.initial.wfn").read_bytes() == previous_files[permanent / "N_3s.initial.wfn"]
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
            ("N_desc.parquet", "rcsf_desc.parquet"),
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

    # Request roots in two blocks with equal J and opposite parity.
    even = [Record((Shell(label, 2),), 0) for label in ("1s", "2s")]
    odd = [
        Record((Shell(label, 1, 1, 1), Shell("2p-", 1, 1, 1)), 0, -1)
        for label in ("1s", "2s")
    ]
    basis = permanent / "blocks.c"
    basis.write_text(
        "Core subshells:\n\nPeel subshells:\n1s 2s 2p-\nCSF(s):\n"
        + "".join(map(text_record, even)) + " *\n" + "".join(map(text_record, odd))
    )
    descriptors(basis)
    basis.unlink()
    run("rwfnestimate", "blocks\ny\n2\n*\n")
    run("rangular" + suffix, "blocks\n" + ("n\n" if launcher else "y\n"))

    def check_roots(path: Path) -> None:
        mixing = records(path)
        assert len(mixing) == 10, path
        assert struct.unpack("=5i", mixing[2]) == (1, 2, 2, 1, 1), path
        assert struct.unpack("=5i", mixing[6]) == (2, 2, 1, 1, -1), path

    root_requirement = "Block 2 (0-) requires at least one valid ASF root"
    missing_root = "Missing ASF selection. " + root_requirement
    blank_selections = ["\n", " \t\n"]
    for solver in ("rmcdhf", "rmcdhf_mem"):
        run(solver + suffix, "blocks\ny\n1-2\n1\n1\n1s\n\n1\n")
        path = permanent / "blocks.mix"
        check_roots(path)
        original = path.read_bytes()
        for blank_selection in blank_selections:
            run(solver + suffix, "blocks\ny\n1\n" + blank_selection, failure=missing_root)
            assert path.read_bytes() == original

    ci_answers = "y\nblocks\n" + ("n\n" if launcher else "") + "n\n" * 6
    run("rci" + suffix, ci_answers + "1\t2\n1\n")
    path = permanent / "blocks.cimix"
    check_roots(path)
    outputs = {item: item.read_bytes() for item in [path, *root.rglob("rci*.res")]}
    assert len(outputs) == (3 if launcher else 2)
    run("rci" + suffix, ci_answers + "\n1\n", failure="Missing ASF selection. Block 1 (0+)")
    assert all(item.read_bytes() == previous for item, previous in outputs.items())
    bad_choices = [(choice, missing_root) for choice in blank_selections]
    bad_choices += [(choice, root_requirement) for choice in ("0\n", "1-3\n", "1-\n")]
    if not launcher:
        # Open MPI's stdin forwarding may leave rank 0 waiting after the launcher sees EOF.
        bad_choices.append(("", missing_root))
    for bad_choice, expected_error in bad_choices:
        run("rci" + suffix, ci_answers + "1\n" + bad_choice, failure=expected_error)
        assert all(item.read_bytes() == previous for item, previous in outputs.items())
    run("rci" + suffix, "n\nblocks\ny\nn\n1-2\n1\n")
    check_roots(path)

    # Mark the disposable MCP headers as diagonal-only to exercise the legacy EAL rejection.
    if launcher:
        for shard in scratch.glob("*/mcp.mcpidx"):
            data = bytearray(shard.read_bytes())
            flags_offset = struct.calcsize("=8s8i8q")
            struct.pack_into("=q", data, flags_offset, 1)
            shard.write_bytes(data)
    else:
        for channel in permanent.glob("mcp.[0-9][0-9]"):
            data = records(channel)
            flags_index = 5 if channel.name == "mcp.30" else 2
            data[flags_index] = struct.pack("=i", 1) + data[flags_index][4:]
            write_records(channel, data)
    for solver in ("rmcdhf", "rmcdhf_mem"):
        run(solver + suffix, "blocks\ny\n1\n1\n1\n1s\n\n1\n",
            failure="Regenerate diagonal-only MCP")


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
