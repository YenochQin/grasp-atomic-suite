!**************************************************************************
!
      subroutine dmerge(n, db, dc, idy, da, dconst, dl)
!-----------------------------------------------
!
!  this merge version has the advantage of loading da(i)
!  and idy(i) only once.
!
!**************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)    :: n
      real(real64), intent(in)    :: dconst
      real(real64), intent(out)   :: dl
      integer, intent(in) :: idy(n)
      real(real64), intent(in) :: db(*)
      real(real64), intent(inout) :: dc(*)
      real(real64), intent(in) :: da(n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: dsum
!-----------------------------------------------

      dsum = 0.0
      do i = 1, n
         dsum = dsum + da(i)*db(idy(i))
         dc(idy(i)) = dc(idy(i)) + dconst*da(i)
      end do
      dl = dsum

      return
      end subroutine dmerge
