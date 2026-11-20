"""Exercise the installed-layout command on disposable RMCDHF run directories."""

import subprocess
import sys
import tempfile
from pathlib import Path


def check_archive(script: Path, run: Path) -> None:
    io_dir = run / "MPI permanent files"
    io_dir.mkdir()
    inputs = {
        "rcsf_desc.parquet": ("case_desc.parquet", b"PAR1\0\xffPAR1"),
        "rcsf_header.toml": ("case_header.toml", b"# byte-bound header\r\n"),
        "rwfn.out": ("case.wfn", b"G92RWF\0"),
        "rmix.out": ("case.mix", b"G92MIX\0"),
        "rmcdhf.sum": ("case.rmcdhf.summary", b" EOL calculation.\n"),
        "isodata": ("isodata", b"nucleus\n"),
        "npol.inp": ("npol.inp", b"nuclear polarization\n"),
    }
    for source, (_, contents) in inputs.items():
        (io_dir / source).write_bytes(contents)
    (run / "rmcdhf.log").write_text("solver log\n")
    (run / "rangular.log").write_text("angular log\n")
    (io_dir / "rcsf.inp").write_text("unused legacy CSFs\n")
    (io_dir / "rcsf_descriptors.parquet").write_bytes(b"older descriptor")

    archive = run / "saved states" / "case"
    command = [script, "--io-dir", io_dir, archive]
    subprocess.run(command, cwd=run, check=True)
    for source, (saved, contents) in inputs.items():
        assert (
            (archive / saved).read_bytes() == (io_dir / source).read_bytes() == contents
        )
    for program in ("rmcdhf", "rangular"):
        assert (archive / f"case.{program}.log").read_bytes() == (
            run / f"{program}.log"
        ).read_bytes()
    assert (archive / "case.w").readlink() == Path(inputs["rwfn.out"][0])
    assert (archive / "case.m").readlink() == Path(inputs["rmix.out"][0])
    assert not (archive / "case.c").exists()

    (io_dir / "rwfn.out").write_bytes(b"next calculation")
    assert subprocess.run(command, cwd=run, check=False).returncode != 0
    assert (archive / "case.w").read_bytes() == inputs["rwfn.out"][1]

    # Serial layout. An average-level run can leave an earlier run's rmix.out in place.
    (io_dir / "rmcdhf.sum").write_text(" (E)AL calculation.\n")
    (io_dir / "rmcdhf.log").touch()
    (io_dir / "npol.inp").unlink()
    (io_dir / "rcsf_desc.parquet").replace(io_dir / "rcsf_descriptors.parquet")
    subprocess.run([script, "average"], cwd=io_dir, check=True)
    average = io_dir / "average"
    assert (average / "average.w").is_file()
    assert not (average / "average.m").is_symlink()
    assert not (average / "average.mix").exists()
    assert not (average / "average.rangular.log").exists()
    assert (average / "average_desc.parquet").read_bytes() == inputs["rcsf_desc.parquet"][1]
    (io_dir / "rcsf_descriptors.parquet").rename(io_dir / "rcsf_desc.parquet")

    # An OL archive needs mixing. A legacy text file cannot replace half of a descriptor pair.
    (io_dir / "rmcdhf.sum").write_text(" OL calculation.\n")
    partial = io_dir / "partial"
    for source in ("rmix.out", "rcsf_desc.parquet", "rcsf_header.toml"):
        (io_dir / source).unlink()
        result = subprocess.run(
            [script, partial], cwd=io_dir, capture_output=True, text=True, check=False
        )
        assert result.returncode != 0 and source in result.stderr, result
        assert not partial.exists()
        assert not list(io_dir.glob(".save_rmcdhf.*"))
        (io_dir / source).write_bytes(inputs[source][1])
    subprocess.run([script, partial], cwd=io_dir, check=True)
    assert (partial / "partial.m").read_bytes() == inputs["rmix.out"][1]


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as directory:
        check_archive(Path(sys.argv[1]).resolve(), Path(directory))
