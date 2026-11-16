      module edensityfit_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine edensityfit (xvec,yvec,z,par,nrnuc,f,rho,res)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points, radial_workspace_capacity
      implicit none (type, external)
      real(real64) :: xvec(radial_workspace_capacity)
      real(real64) :: yvec(max_radial_points)
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: XVEC, YVEC
      real(real64) :: p(5), f(5)
      real(real64) :: par(2)
      real(real64) :: z, drms, rho, res
      integer      :: nrnuc, nparfit
      end subroutine edensityfit
      end interface
      end module edensityfit_i
