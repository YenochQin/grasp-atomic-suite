!**************************************************************************
!
      subroutine dmerge(n, db, dc, idy, da, dconst, dl)
!-----------------------------------------------
!
!  this merge version has the advantage of loading da(i)
!  and idy(i) only once.
!
!**************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:20   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)    :: n
      real(double), intent(in)    :: dconst
      real(double), intent(out)   :: dl
      integer, dimension(n), intent(in) :: idy
      real(double), dimension(*), intent(in) :: db
      real(double), dimension(*), intent(inout) :: dc
      real(double), dimension(n), intent(in) :: da
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: dsum
!-----------------------------------------------

      dsum = 0.0
      do i = 1, n
         dsum = dsum + da(i)*db(idy(i))
         dc(idy(i)) = dc(idy(i)) + dconst*da(i)
      end do
      dl = dsum

      return
      end subroutine dmerge
