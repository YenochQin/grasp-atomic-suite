      module iniest_i
      interface
!
      subroutine iniest(n, nb, niv, hmx, jcol, irow, basis, iblock, jblock)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,only: double
      integer, intent(in) :: n
      integer, intent(in) :: nb
      integer  :: niv
      integer, intent(in) :: jblock
      integer, intent(in) :: jcol(0:*)
      integer, intent(in) :: irow(*)
      integer, intent(in) :: iblock(*)
      real(double), intent(in) :: hmx(*)
      real(double)  :: basis(*)
      end subroutine
      end interface
      end module
