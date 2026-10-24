"""Disposable software fixtures only; never invoke a full atomic solver."""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


SOURCE = "artificial software fixture; E2 convention: test only"
FIELDS = {
    "schema_version": "1", "enabled": ".true.",
    "use_e1": ".true.", "use_e2": ".true.", "units": "'atomic'",
    "model_id": "'flambaum2021_scalar_e1_e2'", "parameter_source": repr(SOURCE),
    "isotope_z": "2", "isotope_a": "4",
    "k1": "16", "b1": "2", "k2": "128", "b2": "2",
}


def config(**changes: str | None) -> str:
    fields = FIELDS | changes
    return "&npol\n" + ",\n".join(f"{key}={value}" for key, value in fields.items() if value is not None) + "\n/\n"


def run(command: list[str], cwd: Path, *, failure: str | None = None) -> str:
    env = os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    result = subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True, timeout=20)
    output = result.stdout + result.stderr
    if failure is None:
        assert result.returncode == 0, (command, result.returncode, output)
        assert "PASS" in output, output
    else:
        assert result.returncode != 0 and failure in output, (command, result.returncode, failure, output)
    return output


def inputs(binary: str, root: Path) -> None:
    path = root / "npol.inp"
    run([binary, "input", str(path)], root)
    valid = [
        "&npol schema_version=1 /", "&npol schema_version=1, enabled=.false., b1=-2 /",
        config(), config(k1="0", k2="0"),
        config(use_e2=".false.", k2=None, b2=None),
        config(use_e1=".false.", k1=None, b1=None),
        config(parameter_source="'it''s an artificial ! fixture / with escaped quote'"),
        "! leading comment\n&NPOL\t schema_version=1 / ! trailing comment\n",
    ]
    for text in valid:
        path.write_text(text)
        run([binary, "input", str(path)], root)
    invalid: list[tuple[str, str]] = [
        ("", "input size"), ("&npol enabled=.false. /", "schema_version"),
        ("&npol schema_version=2 /", "schema_version"),
        ("&npol schema_version=1, mystery=3 /", "invalid namelist"),
        ("&npol schema_version=1", "invalid namelist"),
        ("garbage\n&npol schema_version=1 /", "invalid namelist"),
        ("&npol schema_version=1 / unknown=3", "invalid namelist"),
        ("&npol schema_version=1 /\n&npol schema_version=2 /", "invalid namelist"),
        (config(k1="'bad'"), "invalid namelist"),
        (config(isotope_z="2.5"), "invalid namelist"),
        (config(use_e1=".false.", use_e2=".false."), "use_e1 or use_e2"),
        (config(units=None), "units"), (config(units="'fm'"), "units"),
        (config(model_id="'unknown'"), "model_id"),
        (config(parameter_source=None), "parameter_source"),
        (config(isotope_z=None), "isotope_z"), (config(isotope_z="3"), "isotope_z"),
        (config(isotope_a=None), "isotope_a"), (config(isotope_a="1"), "isotope_a"),
        (config(units=repr("x" * 17)), "units exceeds"),
        (config(model_id=repr("x" * 65)), "model_id exceeds"),
        (config(parameter_source=repr("x" * 1025)), "parameter_source exceeds"),
    ]
    for field in ("k1", "k2", "b1", "b2"):
        for value in (None, "-1", "NaN", "Inf", "-Inf"):
            invalid.append((config(**{field: value}), field))
    for field in ("b1", "b2"):
        invalid.append((config(**{field: "0"}), field))
    for text, error in invalid:
        path.write_text(text)
        run([binary, "input", str(path)], root, failure=error)
    path.write_text(config() + "!" + "x" * 65536)
    run([binary, "input", str(path)], root, failure="input size")
    path.write_text(config())
    path.chmod(0)
    try:
        if os.geteuid() != 0:
            run([binary, "input", str(path)], root, failure="cannot open input")
    finally:
        path.chmod(0o600)
    path.unlink()
    path.mkdir()
    run([binary, "input", str(path)], root, failure="npol:")


def output(binary: str, root: Path) -> None:
    path = root / "npol.inp"
    summary = root / "npol.summary"
    for text in (None, "&npol schema_version=1, enabled=.false. /"):
        if text is not None:
            path.write_text(text)
        run([binary, "init", str(path)], root)
        assert not summary.exists()
    path.write_text(config())
    run([binary, "init", str(path)], root)
    text = summary.read_text()
    for value in ("schema_version = 1", "enabled = T", "use_e1 = T", "use_e2 = T", "isotope_z = 2",
                  "isotope_a = 4", "units = atomic", "model_id = flambaum2021_scalar_e1_e2",
                  f"parameter_source = {SOURCE}", "N = 3", "composed = T"):
        assert value in text, (value, text)
    for key, expected in {"k1 [Eh a0^4]": 16, "b1 [a0]": 2, "k2 [Eh a0^6]": 128, "b2 [a0]": 2}.items():
        line = next(line for line in text.splitlines() if line.startswith(key))
        assert float(line.split("=", 1)[1]) == expected
    summary.unlink()
    path.write_text(config(k1="0", k2="0"))
    run([binary, "init", str(path)], root)
    summary.unlink()
    # Replacing a longer previous summary cannot leave old metadata at EOF.
    summary.write_text("old summary" * 2000)
    path.write_text(config(use_e2=".false.", k2=None, b2=None))
    run([binary, "init", str(path)], root)
    assert "use_e2 = F" in summary.read_text()
    assert summary.read_text().endswith("composed = T\n")
    summary.unlink()
    summary.mkdir()
    run([binary, "init", str(path)], root, failure="cannot open summary")
    summary.rmdir()
    if Path("/dev/full").exists():
        summary.symlink_to("/dev/full")
        run([binary, "init", str(path)], root, failure="cannot write/close summary")
        assert not summary.exists()
    path.write_text(config(k1="NaN"))
    run([binary, "init", str(path)], root, failure="k1")
    assert not summary.exists()
    path.write_text(config(b1="1e-300"))
    run([binary, "init", str(path)], root, failure="not representable")
    assert not summary.exists()


def entrypoints(root: Path) -> None:
    for variant in ("rmcdhf", "rmcdhf_mem", "rmcdhf_mpi", "rmcdhf_mem_mpi"):
        directory = root / "src" / "appl" / variant
        mpi = variant.endswith("_mpi")
        stem = "getscdmpi" if mpi else "getscd"
        hook = "npol_init_mpi" if mpi else "npol_init"
        text = (directory / f"{stem}.f90").read_text().lower()
        # Ignore comments when checking code ordering and call counts.
        code = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
        calls = re.findall(r"call\s+(\w+)", code)
        assert calls.count(hook) == 1
        assert calls.index("nucpot") < calls.index(hook) < calls.index("setrwfmpi" if mpi else "setrwfa")
        for filename in (f"{stem}.f90", f"{stem}_I.f90"):
            declaration = (directory / filename).read_text().lower()
            assert re.search(r"subroutine\s+" + stem + r"\s*\(eol, idblk, isofile, rwffile, npfile\)", declaration)
        main = (directory / ("rscfmpivu.f90" if mpi else "rscfvu.f90")).read_text().lower()
        if mpi:
            assert "npfile = ''" in main and "npfile = trim(permdir) // '/npol.inp'" in main
            assert re.search(r"call getscdmpi\s*\(eol, idblk, file1, file2, npfile\)", main)
        else:
            assert "call getscd (eol, idblk, 'isodata', 'rwfn.inp', 'npol.inp')" in main
        # No second matrix contribution and MODE=0 retained in both layouts.
        setham = (directory / "setham.f90").read_text().lower()
        assert "npol" not in setham and re.search(r"rinti\s*\([^\n]*,\s*0\s*\)", setham)
    assert "npol" not in (root / "src/lib/lib9290/nucpot.f90").read_text().lower()
    for variant in ("rmcdhf90", "rmcdhf90_mem", "rmcdhf90_mpi", "rmcdhf90_mem_mpi",
                    "rmcdhf90_mem_strict_mpi", "rmcdhf90_mem_diag_mpi"):
        stem = "getscdmpi" if variant.endswith("_mpi") else "getscd"
        text = (root / "src/appl" / variant / f"{stem}.f90").read_text().lower()
        assert "npol_init" not in text and "npfile" not in text, variant


def mpi_cases(command: list[str], root: Path, failure: bool) -> None:
    scenarios = ("parse", "length", "mismatch", "charge", "compose", "summary", "build") if failure else (
        "ok", "disabled", "missing", "zero")
    for mode in scenarios:
        case = root / mode
        for child in ("rank0", "rank1", "input"):
            (case / child).mkdir(parents=True)
        (case / "rank1/worker-must-not-read.inp").write_text("INVALID WORKER FILE")
        path = case / "input/npol.inp"
        if mode != "missing":
            path.write_text({"parse": "&npol schema_version=99 /", "disabled": "&npol schema_version=1 /",
                             "zero": config(k1="0", k2="0"), "build": config(b1="1e-300")}.get(mode, config()))
        if mode == "summary":
            (case / "input/npol.summary").mkdir()
        run(command + [mode, str(case)], case, failure="EXPECTED_COLLECTIVE_FAILURE" if failure else None)
        assert not (case / "rank0/npol.summary").exists()
        assert not (case / "rank1/npol.summary").exists()
        summary = case / "input/npol.summary"
        if mode in ("ok", "zero"):
            assert "composed = T" in summary.read_text()
        elif mode != "summary":
            assert not summary.exists()


def main() -> None:
    kind, *command = sys.argv[1:]
    with tempfile.TemporaryDirectory(prefix="npol-code-") as directory:
        root = Path(directory)
        if kind == "input":
            inputs(command[0], root)
        elif kind == "output":
            output(command[0], root)
        elif kind == "entrypoints":
            entrypoints(Path(command[0]))
            run([command[1], "compose"], root)
        elif kind in ("mpi", "mpi_error"):
            mpi_cases(command, root, kind == "mpi_error")
        elif kind == "run":
            run(command, root)
        else:
            raise ValueError(kind)
    print(f"PASS {kind}")


if __name__ == "__main__":
    main()
