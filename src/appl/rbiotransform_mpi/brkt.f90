!***********************************************************************
!                                                                      *
      subroutine brkt
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
      use biorb_c
      use wave_c, only:  mfii, pfii, qfii, mfff, pfff, qfff
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, l
      real(double), dimension(nnnw,nnnw) :: braket
      real(double) :: result
!-----------------------------------------------
!

      do i = 1, nwii
         do j = 1, nwff
            if (nakii(i) /= nakff(j)) cycle
!
!   Determine the maximum tabulation point for the integrand
!
            mtp = min(mfii(i),mfff(j))
!
!   Tabulate the integrand as required for SUBROUTINE QUAD; the
!   value at the first tabulation point is arbitrary
!
            ta = 0.0d0
!            TA(1)=0.D0
            do l = 2, mtp
               ta(l) = (pfii(l,i)*pfff(l,j) + qfii(l,i)*qfff(l,j))*rp(l)
            end do
!
!   Perform the quadrature
!
            call quad (result)

            braket(i,j) = result
            if (myid .eq. 0)                                          &
            write (*,9) '<',npii(i),nhii(i),'|',npff(j),nhff(j),'> =',&
               braket(i,j)
         end do
      end do
!
    9 format(a,i2,a,a,i2,a,a,e20.13)

      return
      end subroutine brkt
