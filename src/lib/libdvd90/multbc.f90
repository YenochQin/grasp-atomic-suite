module multbc_i
    implicit none (type, external)
    private
    public :: multbc
contains
      subroutine multbc(n, k, m, c, temp, b)
!=======================================================================
!     called by: DVDRVR
!
!     Multiplies B(N,K)*C(K,M) and stores it in B(N,M)
!     Used for collapsing the expanding basis to current estimates,
!     when basis becomes too large, or for returning the results back
!
!       Subroutines called
!       DINIT, DGEMV, DCOPY
!-----------------------------------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dgemv_I
      !USE dcopy_I
      use suite_linear_algebra, only: dcopy, dgemv
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: k
      integer, intent(in) :: m
      real(real64)  :: c(k*m)
      real(real64)  :: temp(m)
      real(real64)  :: b(n*k)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: irow
!-----------------------------------------------
!-----------------------------------------------------------------------
      do irow = 1, n
         call dgemv ("Transp", k, m, 1.0_real64, c, k, b(irow), n, 0.0_real64, temp, 1)
         call dcopy (m, temp, 1, b(irow), n)
      end do

      return
      end subroutine multbc
end module multbc_i
