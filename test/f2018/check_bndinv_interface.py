"""Check both BNDINV module procedures for explicit bounds and the public calling contract."""
from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys
import tempfile


PROBE = """program bndinv_contract
    use, intrinsic :: iso_fortran_env, only: real64
    use invert_banded_matrix_interface, only: invert_banded_matrix
    implicit none (type, external)
    abstract interface
        subroutine banded_inverse(matrix, inverse_workspace, matrix_order, determinant, singularity_tolerance, &
            inversion_status, leading_dimension)
            import :: real64
            integer, intent(in) :: matrix_order, leading_dimension
            real(real64), intent(inout) :: matrix(leading_dimension, leading_dimension)
            real(real64), intent(inout) :: inverse_workspace(leading_dimension, leading_dimension)
            real(real64), intent(out) :: determinant
            real(real64), intent(in) :: singularity_tolerance
            integer, intent(out) :: inversion_status
        end subroutine banded_inverse
    end interface
    procedure(banded_inverse), pointer :: inverse
    real(real64) :: matrix(1, 1), inverse_workspace(1, 1), determinant
    integer :: inversion_status

    inverse => invert_banded_matrix
    call invert_banded_matrix(leading_dimension=1, inversion_status=inversion_status, singularity_tolerance=0.0_real64, &
        determinant=determinant, matrix_order=1, inverse_workspace=inverse_workspace, matrix=matrix)
end program bndinv_contract
"""


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    for variant in ("biorthtransform", "biorthtransform_mpi"):
        app = root / "src/apps" / variant
        implementation = app / "invert_banded_matrix.f90"
        source = implementation.read_text()
        # GNU's procedure association check does not compare explicit bounds.
        # Check these separately instead of treating compilation as proof.
        for name in ("matrix", "inverse_workspace"):
            pattern = rf"(?:^|::|,)\s*{name}\s*\(([^)]+)\)"
            actual = re.search(pattern, source, re.I | re.M)
            assert actual is not None and re.sub(r"\s+", "", actual.group(1)).lower() == "leading_dimension,leading_dimension", (variant, name)
        with tempfile.TemporaryDirectory(prefix="bndinv-interface-") as temporary:
            work = Path(temporary)
            (work / "probe.f90").write_text(PROBE)
            subprocess.run([sys.argv[1], "-std=f2018", "-pedantic-errors", "-fimplicit-none",
                            "-fsyntax-only", str(implementation), "probe.f90"],
                           cwd=work, check=True, capture_output=True, text=True)
        print(f"{variant}: BNDINV bounds, procedure association and keyword arguments passed")


if __name__ == "__main__":
    main()
