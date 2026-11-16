      module iniestmpi_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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
