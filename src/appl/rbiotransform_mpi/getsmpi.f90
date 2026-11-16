!***********************************************************************
!                                                                      *
      subroutine gets(s, nshli, nshlf)
!                                                                      *
!   This subroutine calculates the initial and final state             *
!   radial overlap matrix                                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use grid_c
      use tatb_c
      use orb_c
      use biorb_c
      use wave_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      use wrtmat_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer       :: nshli
      integer       :: nshlf
      real(double)  :: s(nshli,nshlf)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: i, j, l
      real(double) :: result
!-----------------------------------------------
!
      do i = 1, nshli
         do j = 1, nshlf
!
!   Determine the maximum tabulation point for the integrand
!
            mtp = min(mfii(i),mfff(j))
!
!   Tabulate the integrand as required for SUBROUTINE QUAD; the
!   value at the first tabulation point is arbitrary
!
            ta(1) = 0.d0
            do l = 2, mtp
               ta(l) = (pfii(l,i)*pfff(l,j) + qfii(l,i)*qfff(l,j))*rp(l)
            end do
!
!   Perform the quadrature
!
            call quad (result)
            s(i,j) = result
         end do
      end do
!
!   Print out
!
      if (myid .eq. 0) then
!      WRITE (*, *) '********************'
!      WRITE (*, *) ' S matrix from GETS'
!      WRITE (*, *) '********************'

         call wrtmat (s, nshli, nshlf, nshli, nshlf)
      endif
!
      return
      end subroutine gets
