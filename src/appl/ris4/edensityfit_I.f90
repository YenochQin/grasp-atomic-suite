      module edensityfit_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine edensityfit (xvec,yvec,z,par,nrnuc,f,rho,res)
      use vast_kind_param,  only: double
      use parameter_def,    only: max_radial_points, radial_workspace_capacity
      real(double), dimension(radial_workspace_capacity) :: xvec
      real(double), dimension(max_radial_points) :: yvec
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: XVEC, YVEC
      real(double), dimension(5) :: p, f
      real(double), dimension(2) :: par
      real(double) :: z, drms, rho, res
      integer      :: nrnuc, nparfit
      end subroutine
      end interface
      end module
