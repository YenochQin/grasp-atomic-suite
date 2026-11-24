"""Check import_rsave with archived upstream GRASP files."""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def scan(directory: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in directory.iterdir()}


def check_import(script: Path, gj: Path, data: Path, rcsfs: str, root: Path) -> None:
    command = [sys.executable, script, "--rcsfs", rcsfs, "test"]
    source_to_target = {
        "test.w": "test.wfn",
        "test.m": "test.mix",
        "test.sum": "test.rmcdhf.summary",
        "test.log": "test.rmcdhf.log",
        "test.alog": "test.rangular.log",
    }
    for name in (*source_to_target, "test.c", "isodata", "test.lsj.lbl"):
        shutil.copy2(data / name, root / name)
    original = scan(root)
    subprocess.run(command, cwd=root, check=True, timeout=45)
    imported = scan(root)
    assert all(imported[name] == contents for name, contents in original.items())
    assert all(imported[target] == original[source] for source, target in source_to_target.items())
    assert set(imported) - set(original) == {*source_to_target.values(), "test_desc.parquet", "test_header.toml"}

    result = subprocess.run(command, cwd=root, check=False, capture_output=True, text=True, timeout=45)
    assert result.returncode != 0 and "destination already exists: test.wfn" in result.stderr, result
    assert scan(root) == imported

    result = subprocess.run(
        [gj], cwd=root, input="y\ntest\nn\n", check=True, capture_output=True, text=True, timeout=45
    )
    assert "Execution complete" in result.stdout + result.stderr, result
    assert (root / "test.gj").read_bytes() == (data / "test.gj").read_bytes()

    average = root / "average"
    average.mkdir()
    for name in ("test.c", "test.w", "test.m"):
        shutil.copyfile(data / name, average / name)
    (average / "test.sum").write_text(" (E)AL calculation.\n")
    # A separately saved CI file uses the same binary transport, with its own suffix.
    (average / "test.cm").write_bytes(b"separate CI mixing records")
    subprocess.run(command, cwd=average, check=True, timeout=45)
    assert not (average / "test.mix").exists()
    assert (average / "test.m").read_bytes() == (data / "test.m").read_bytes()
    assert (average / "test.cimix").read_bytes() == (average / "test.cm").read_bytes()

    invalid = root / "invalid"
    invalid.mkdir()
    (invalid / "test.c").write_text("invalid CSFs\n")
    shutil.copyfile(data / "test.w", invalid / "test.w")
    before = scan(invalid)
    result = subprocess.run(command, cwd=invalid, check=False, capture_output=True, text=True, timeout=45)
    assert result.returncode != 0 and "rcsfs convert failed" in result.stderr, result
    assert scan(invalid) == before


if __name__ == "__main__":
    executable = os.environ.get("RCSFS_TEST_CLI") or shutil.which("rcsfs")
    if executable is None:
        print("V3-capable rcsfs is required by this fixture")
        raise SystemExit(77)
    script, gj, data = (Path(value).resolve() for value in sys.argv[1:])
    with tempfile.TemporaryDirectory(prefix="import-rsave-") as directory:
        check_import(script, gj, data, executable, Path(directory))
