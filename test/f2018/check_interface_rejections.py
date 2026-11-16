"""Prove the public interfaces accept valid calls and reject real ABI errors."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> None:
    compiler, moddir, mpidir = sys.argv[1:4]
    includes = ["-I" + moddir, "-I" + mpidir, *sys.argv[4:]]
    common = """program probe
 use, intrinsic :: iso_fortran_env, only: real64, int64
 use suite_linear_algebra
 implicit none (type, external)
 real(real64), pointer :: p(:)
 real(real64) :: y(3)
 {declaration}
 allocate(p(3))
 p=1
 {call}
end program
"""
    cases = {
        "valid_vector": (common.format(declaration="", call="call dcopy(3,p(:),1,y,1)"), True),
        "pointer_element": (common.format(declaration="", call="call dcopy(3,p(1),1,y,1)"), False),
        "wrong_count_kind": (common.format(declaration="", call="call dcopy(3_int64,p(:),1,y,1)"), False),
        "wrong_return_type": (common.format(declaration="integer, external :: ddot", call="y(1)=ddot(3,p(:),1,p(:),1)"), False),
    }
    for kind, valid in [("integer", True), ("integer(int8)", False)]:
        cases["mpi_handle_" + str(valid)] = (
            "program probe\n use mpi_C\n use, intrinsic :: iso_fortran_env, only: int8\n implicit none (type, external)\n"
            + kind + " :: dtype\n call mpix_bytes(1,dtype,mpi_error_code)\n end program\n", valid
        )
    with tempfile.TemporaryDirectory(prefix="suite-interface-") as temporary:
        work = Path(temporary)
        for name, (source, valid) in cases.items():
            path = work / f"{name}.f90"
            path.write_text(source)
            result = subprocess.run([compiler, "-std=f2018", "-fsyntax-only", *includes, str(path)],
                                    cwd=work, capture_output=True, text=True, timeout=30)
            assert (result.returncode == 0) == valid, (name, result.stdout, result.stderr)
            print(f"{name}: {'accepted' if valid else 'rejected'}")


if __name__ == "__main__":
    main()
