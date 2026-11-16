!***********************************************************************
      subroutine orbout(rwffile2)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use grid_c
      use orb_c
      use wave_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: rwffile2*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, mfj, i
!-----------------------------------------------
!
      open(23, file=rwffile2, status='UNKNOWN', form='UNFORMATTED', position=&
         'asis')
      write (23) 'G92RWF'
      do j = 1, nw
         mfj = mf(j)
         write (23) np(j), nak(j), e(j), mfj
         write (23) pz(j), (pf(i,j),i=1,mfj), (qf(i,j),i=1,mfj)
         write (23) (r(i),i=1,mfj)               ! This is a waste of resources
      end do
      close(23)

      return
      end subroutine orbout
