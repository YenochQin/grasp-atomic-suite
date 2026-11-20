"""Run both complete RATDEN bodies with controlled orbital/quadrature seams."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> None:
    compiler, modules, library = sys.argv[1:4]
    root = Path(sys.argv[4]) if len(sys.argv) > 4 else Path(__file__).resolve().parents[2]
    support = """module dirac_coulomb_radial_orbital
use iso_fortran_env, only: real64
use grasp_parameters, only: max_radial_points
implicit none
contains
subroutine dcbsrw(n,kappa,z,e,rg0,rg,rf,mtp)
integer,intent(in)::n,kappa
integer,intent(out)::mtp
real(real64),intent(in)::z
real(real64),intent(out)::e,rg0,rf(max_radial_points)
real(real64),intent(inout)::rg(max_radial_points)
e=0; rg0=0; rg=1; rf=0; mtp=4
end subroutine
end module
module radial_quadrature
use iso_fortran_env, only: real64
use radial_quadrature_workspace, only: mtp,ta
implicit none
contains
subroutine quad(result)
real(real64),intent(out)::result
result=sum(ta(1:mtp))
end subroutine
end module
"""
    driver = """program qed_boundary
use iso_fortran_env, only: real64,real32
use grasp_parameters, only: max_radial_points
use radial_grid_state, only: radial_grid_derivative
use radial_quadrature_workspace, only: mtp
use ratden_probe
implicit none
real(real64)::p(max_radial_points),q(max_radial_points),z,value,points(5)
integer::i,n,k
n=1; k=-1; z=1; p=1; q=0
points=[nearest(0.0219_real64,-1.0_real64),0.0219_real64, &
nearest(0.0219_real64,1.0_real64), &
(real(0.0219_real32,real64)+0.0219_real64)/2,0.01_real64]
do i=1,size(points)
radial_grid_derivative=1
radial_grid_derivative(2)=0.001_real64
radial_grid_derivative(3)=points(i)
value=ratden(p,q,4,n,k,z)
if(mtp /= merge(3,2,points(i)<=0.0219_real64)) error stop 'RATDEN cutoff boundary'
if(value/=1.0_real64) error stop 'RATDEN controlled ratio'
enddo
end program
"""
    for variant in ("rci", "rci_mpi"):
        text = (root / f"src/appl/{variant}/qed.f90").read_text()
        start = text.index("      real(real64) function ratden")
        end = text.index("      end function ratden", start) + len("      end function ratden")
        with tempfile.TemporaryDirectory(prefix="qed-boundary-") as temporary:
            work = Path(temporary)
            (work / "support.f90").write_text(support)
            (work / "ratden.f90").write_text("module ratden_probe\ncontains\n" + text[start:end] + "\nend module\n")
            (work / "driver.f90").write_text(driver)
            subprocess.run([compiler, "-std=f2018", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
                            "-I" + modules, "support.f90", "ratden.f90", "driver.f90", library, "-o", "probe"],
                           cwd=work, check=True, capture_output=True, text=True)
            subprocess.run([work / "probe"], cwd=work, check=True)
        print(f"{variant}: RATDEN exact, adjacent, midpoint and normal grid cutoff passed")


if __name__ == "__main__":
    main()
