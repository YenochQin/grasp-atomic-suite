"""Compile the actual EPS declarations and private environment helper in isolation."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("compiler")
    parser.add_argument("module_directory", type=Path)
    parser.add_argument("--baseline-ref")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]

    def source(relative: str) -> str:
        if args.baseline_ref:
            return subprocess.check_output(["git", "show", f"{args.baseline_ref}:{relative}"], cwd=root, text=True)
        return (root/relative).read_text()

    with tempfile.TemporaryDirectory(prefix="input-contract-") as temporary:
        work = Path(temporary)
        for variant in ("rmcdhf", "rmcdhf_mem", "rmcdhf_mpi", "rmcdhf_mem_mpi"):
            text = source(f"src/appl/{variant}/setcof.f90")
            statements = re.findall(r"^\s*real[^\n]*,\s*parameter\s*::\s*eps\s*=[^\n]+", text, re.I | re.M)
            assert len(statements) == 1
            probe = ("program eps_contract\nuse vast_kind_param, only: double\nimplicit none\n"
                     + statements[0] + "\nif (kind(eps) /= double) error stop 'EPS is not DOUBLE'\n"
                     + "if (real(eps,double) /= 1.0d-10) error stop 'EPS value lost precision'\nend program\n")
            path = work/f"{variant}.f90"
            path.write_text(probe)
            binary = work/variant
            subprocess.run([args.compiler, "-std=f2018", "-fcheck=all", "-I"+str(args.module_directory.resolve()),
                            str(path), "-o", str(binary)], check=True, cwd=work, capture_output=True, text=True)
            result = subprocess.run([str(binary)], capture_output=True, text=True)
            assert result.returncode == 0, (variant, result.stdout+result.stderr)
            print(variant + ": DOUBLE EPS verified")

        text = source("src/lib/mpi90/cpath.f90")
        helper = re.search(r"^\s*function getenv_allocating\([^\n]+\n.*?^\s*end\s*$", text, re.M | re.S)
        assert helper
        # The contained helper is private to CPATH; compile its exact source
        # body, without changing production visibility or introducing a copy.
        program = """program environment_contract
implicit none
character(:), allocatable :: value
logical :: found
found = getenv_allocating('GRASP_IO_TEST_MISSING', value)
if (found .or. allocated(value)) error stop 'missing variable'
found = getenv_allocating('GRASP_IO_TEST_SHORT', value)
if (.not.found) error stop 'first allocation'
if (value /= 'a') error stop 'short value'
found = getenv_allocating('GRASP_IO_TEST_SHORT', value)
if (.not.found .or. value /= 'a') error stop 'same length'
found = getenv_allocating('GRASP_IO_TEST_LONG', value)
if (.not.found .or. value /= 'long path with spaces') error stop 'reallocation'
found = getenv_allocating('GRASP_IO_TEST_EMPTY', value)
if (.not.found .or. len(value) /= 0) error stop 'empty variable'
found = getenv_allocating('GRASP_IO_TEST_MISSING', value)
if (found .or. len(value) /= 0) error stop 'missing variable changed value'
print *, 'PASS CPATH environment helper'
contains
""" + helper[0] + "\nend program\n"
        path = work/"environment.f90"
        path.write_text(program)
        binary = work/"environment"
        subprocess.run([args.compiler, "-std=f2018", "-fcheck=all", "-fno-frontend-optimize", str(path), "-o", str(binary)],
                       check=True, cwd=work, capture_output=True, text=True)
        env = os.environ | {"GRASP_IO_TEST_SHORT": "a", "GRASP_IO_TEST_LONG": "long path with spaces", "GRASP_IO_TEST_EMPTY": ""}
        env.pop("GRASP_IO_TEST_MISSING", None)
        subprocess.run([str(binary)], env=env, check=True)


if __name__ == "__main__":
    main()
