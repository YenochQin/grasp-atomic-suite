      module iniestmpi_i
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine iniestmpi (nmax, ncf, niv, basis, hmx, jcol, irow)
      use vast_kind_param, only:  double
      integer, intent(in) :: nmax
      integer, intent(in) :: ncf
      integer  :: niv
      integer, intent(in) :: jcol(0:*)
      integer, intent(in) :: irow(*)
      real(double)  :: basis(*)
      real(double), intent(in) :: hmx(*)
      end subroutine
      end interface
      end module
