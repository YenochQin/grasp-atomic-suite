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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grid_c
      use tatb_c
      use orb_c
      use biorb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      use wrtmat_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nshli
      integer, intent(in) :: nshlf
      real(real64)  :: s(nshli,nshlf)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: i, j, l
      real(real64) :: result
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
            ta(1) = 0._real64
            do l = 2, mtp
               ta(l) = (pfii(l,i)*pfff(l,j) + qfii(l,i)*qfff(l,j))*radial_grid_derivative(l)
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
      write (*, *) '********************'
      write (*, *) ' S matrix from GETS'
      write (*, *) '********************'

      call wrtmat (s, nshli, nshlf, nshli, nshlf)
!
      return
      end subroutine gets
