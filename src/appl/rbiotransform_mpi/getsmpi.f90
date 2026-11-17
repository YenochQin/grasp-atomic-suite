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
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use tatb_c, only: double, mtp, radial_workspace_capacity, ta, tb
      use orb_c, only: double, max_orbitals
      use biorb_c, only: double, max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use wave_c, only: double, max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use mpi_c, only: double, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i, only: quad
      use wrtmat_i, only: wrtmat
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
            ta(1) = 0.0_real64
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
      if (mpi_rank == 0) then
!      WRITE (*, *) '********************'
!      WRITE (*, *) ' S matrix from GETS'
!      WRITE (*, *) '********************'

         call wrtmat (s, nshli, nshlf, nshli, nshlf)
      end if
!
      return
      end subroutine gets
