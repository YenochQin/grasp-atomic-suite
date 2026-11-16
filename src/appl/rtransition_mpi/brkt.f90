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
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:56:42   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: max_orbitals
      use grid_c
      use tatb_c
      use biorb_c
      use wave_c, only: mfii, pfii, qfii, mfff, pfff, qfff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, l
      real(double), dimension(max_orbitals,max_orbitals) :: braket
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
            ta(1) = 0.d0
            do l = 2, mtp
               ta(l) = (pfii(l,i)*pfff(l,j) + qfii(l,i)*qfff(l,j))*radial_grid_derivative(l)
            end do
!
!   Perform the quadrature
!
            call quad (result)

            braket(i,j) = result

!            WRITE (*,9) '<',NPII(I),NHII(I),'|',NPFF(J),NHFF(J),'> =', &
!               BRAKET(I,J)
         end do
      end do
!
!    9 FORMAT(A,I2,A,A,I2,A,A,E20.13)

      return
      end subroutine brkt
