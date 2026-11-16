      module iniest2_i
      implicit none (type, external)
      interface
      subroutine iniest2 (nmax, ncf, niv, basis, hmx, jcol, irow)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nmax
      integer, intent(in) :: ncf
      integer, intent(in) :: niv
      real(real64) :: basis(*)
      real(real64), intent(in) :: hmx(*)
      integer, intent(in) :: jcol(0:*)
      integer, intent(in) :: irow(*)
!VAST...Calls: ALLOC, DINIT, DSPEVX, DCOPY, DALLOC
      end subroutine
      end interface
      end module
