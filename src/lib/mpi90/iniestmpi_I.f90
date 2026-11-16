      module iniestmpi_i
      implicit none (type, external)
      interface
      subroutine iniestmpi (nmax, ncf, niv, basis, hmx, jcol, irow)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nmax
      integer, intent(in) :: ncf
      integer  :: niv
      integer, intent(in) :: jcol(0:*)
      integer, intent(in) :: irow(*)
      real(real64)  :: basis(*)
      real(real64), intent(in) :: hmx(*)
      end subroutine
      end interface
      end module
