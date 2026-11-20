"""Compile the reviewed production literals and test double threshold boundaries.

This is an expression-level test, not a substitute for SCF convergence tests.
Exact literals, the SLAMCH ABI and memory-report arithmetic are outside its scope.
"""
from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys
import tempfile


def candidates() -> list[tuple[str, str, int]]:
    cases = [(f"src/apps/{name}/{routine}.f90", "0.001", 1)
             for name, routine in (("rmcdhf", "scf"), ("rmcdhf_mem", "scf"),
                                   ("rmcdhf_mpi", "scfmpi"), ("rmcdhf_mem_mpi", "scfmpi"))]
    cases += [("src/apps/rwfnestimate/genrwf.f90", "0.05", 1),
              ("src/apps/rwfnestimate/solvh.f90", "0.1e00", 1),
              ("src/apps/rwfnestimate/solvh.f90", "0.2e00", 2)]
    cases += [(f"src/apps/{name}/qed.f90", "0.0219", 1) for name in ("rci", "rci_mpi")]
    cases += [(f"src/apps/{name}/bndinv.f90", "1.0e-8", 1) for name in ("rbiotransform", "rbiotransform_mpi")]
    cases += [("src/apps/jj2lsj/jj2lsj_code.f90", literal, 1)
              for literal in ("0.01", "1.0001", "0.00001", "1.001", ".998")]
    cases += [(f"src/apps/{name}/printaLS.f90", "0.00000001", 1)
              for name in ("rtransition", "rtransition_mpi", "rtransition_phase")]
    return cases


def main() -> None:
    compiler = sys.argv[1]
    root = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parents[2]
    statements = []
    count = 0
    for relative, literal, occurrences in candidates():
        text = "\n".join(line.split("!")[0] for line in (root / relative).read_text().splitlines())
        # A leading zero is optional Fortran spelling, not a different threshold.
        spelling = (r"0?" if literal.startswith(".") else "") + re.escape(literal)
        matches = re.findall(r"(?<![\w.])" + spelling + r"_real(?:32|64)\b", text)
        assert len(matches) == occurrences, (relative, literal, matches)
        for token in matches:
            count += 1
            statements += [f"actual = {token}", f"expected = {literal}_real64",
                           f"call check('{relative}: {literal}', actual, expected)"]
    source = """program literal_boundaries
use iso_fortran_env, only: real32, real64
implicit none
real(real64) :: actual, expected
""" + "\n".join(statements) + """
contains
subroutine check(label, actual, expected)
character(*), intent(in) :: label
real(real64), intent(in) :: actual, expected
real(real64) :: scales(3), threshold, reference, points(5), x
integer :: i,j
scales = [1.0_real64, 1e-8_real64, 100.0_real64]
do i=1,size(scales)
    threshold = actual*scales(i)
    reference = expected*scales(i)
    points = [nearest(reference,-1.0_real64), reference, nearest(reference,1.0_real64), &
              reference/2.0_real64, reference*2.0_real64]
    do j=1,size(points)
        x = points(j)
        if ((x < threshold) .neqv. (x < reference)) then
            print *, trim(label), ' strict boundary mismatch', x, threshold, reference
            error stop 1
        endif
        if ((x <= threshold) .neqv. (x <= reference)) then
            print *, trim(label), ' inclusive boundary mismatch', x, threshold, reference
            error stop 1
        endif
    enddo
enddo
end subroutine
end program
"""
    with tempfile.TemporaryDirectory(prefix="numeric-literals-") as temporary:
        work = Path(temporary)
        (work / "probe.f90").write_text(source)
        subprocess.run([compiler, "-std=f2018", "-fcheck=all", "probe.f90", "-o", "probe"],
                       cwd=work, check=True, capture_output=True, text=True)
        subprocess.run([work / "probe"], cwd=work, check=True)
    print(f"{count} production literals: adjacent/exact/typical boundaries at three scales passed")


if __name__ == "__main__":
    main()
