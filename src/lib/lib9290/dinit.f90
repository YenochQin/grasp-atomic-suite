module dinit_i
    implicit none (type, external)
    private
    public :: dinit
contains
      subroutine dinit(n, a, x, incx)
!     ==================================================================
!
!     PURPOSE ... INITIALIZES REAL*8           VECTOR TO
!                 A CONSTANT VALUE 'A'
!
!     CREATED ... APR. 14, 1987
!
!     ==================================================================
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: incx
      real(real64), intent(in) :: a
      real(real64), intent(out) :: x(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: xaddr, i
!-----------------------------------------------
!
      if (incx == 1) then
!
!         ----------------------------------
!         ... UNIT INCREMENT (STANDARD CASE)
!         ----------------------------------
!
         x(:n) = a
!
      else
!
!         ----------------------
!         ... NON-UNIT INCREMENT
!         ----------------------
!
         xaddr = 1
         if (incx < 0) xaddr = ((-n) + 1)*incx + 1
!
         x(xaddr:(n-1)*incx+xaddr:incx) = a
!
      end if
!
      return
!
      end subroutine dinit
end module dinit_i
