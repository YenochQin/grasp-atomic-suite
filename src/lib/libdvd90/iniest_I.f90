      module iniest_i
      implicit none (type, external)
      interface
!
      subroutine iniest(n, nb, niv, hmx, jcol, irow, basis, iblock, jblock)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: nb
      integer  :: niv
      integer, intent(in) :: jblock
      integer, intent(in) :: jcol(0:*)
      integer, intent(in) :: irow(*)
      integer, intent(in) :: iblock(*)
      real(real64), intent(in) :: hmx(*)
      real(real64)  :: basis(*)
      end subroutine
      end interface
      end module
