!***********************************************************************
!                                                                      *
      subroutine ltab(is, nqs, ks, irows)
!                                                                      *
!   locates rows of possible parents of active shell states for acc-   *
!   essing  NTAB. It is assumed that empty shells have been elimina-   *
!   ted from consideration by SUBROUTINE RKCO.                         *
!                                                                      *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:49:59   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use terms_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: is(4)
      integer, intent(inout) :: nqs(4)
      integer, intent(in) :: ks(4)
      integer, intent(inout) :: irows(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kq(4)
      integer :: kq1, kq2
      integer  :: i
!-----------------------------------------------
!
!
      if (is(1) == is(2)) nqs(1) = nqs(2) - 1
      if (is(3) == is(4)) nqs(3) = nqs(4) - 1
!
      do i = 1, 4
!
!   Check that input data are consistent
!
         if (nqs(i)<=0 .or. nqs(i)>ks(i)) then
            write (*, 300) nqs(i), is(i), ks(i)
            stop
         endif
!
         kq1 = nqs(i) - 1
         kq2 = ks(i) - kq1
         kq(i) = min(kq1,kq2) + 1
         if (kq(i) /= 1) then
            irows(i) = (ks(i)*(ks(i)-2))/8 + kq(i)
         else
            irows(i) = 1
         endif
!
         if (irows(i) <= nrows) cycle
         write (*, 301)
         stop
!
      end do
!
      return
!
  300 format('LTAB: ',1i3,' Electrons in shell ',1i3,' with 2j+1 = ',1i3)
  301 format('LTAB: Extend COMMON block TERMS')
      return
!
      end subroutine ltab
