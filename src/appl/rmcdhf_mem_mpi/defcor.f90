!***********************************************************************
!                                                                      *
      subroutine defcor(j)
!                                                                      *
!   Compute the deferred corrections for orbital J .                   *
!                                                                      *
!                                          Last updated: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use grid_c
      use wave_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: w3 = 1.0_real64/120.0_real64
      real(real64), parameter :: w2 = -15.0_real64*w3
      real(real64), parameter :: w1 = 40.0_real64*w3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, mfjm3, mfjm2
      logical :: first
!-----------------------------------------------
!
      data first/ .true./
!
!   The deferred corrections for the first two points are
!   unnecessary, because the integration always commences
!   from the fourth point of the grid; this requires the
!   deferred correction at the third and subsequent points
!   only
!
      if (first) then
         dp(:2) = 0.0_real64
         dq(:2) = 0.0_real64
         first = .false.
      endif
!
!   Intermediate points
!
      mfjm3 = orbital_last_grid_index(j) - 3
      do i = 3, mfjm3
!
         dp(i) = w3*(radial_large_components(i + 3,j) - radial_large_components(i - 2,j)) + &
                 w2*(radial_large_components(i + 2,j) - radial_large_components(i - 1,j)) + &
                 w1*(radial_large_components(i + 1,j) - radial_large_components(i,j))
!
         dq(i) = w3*(radial_small_components(i + 3,j) - radial_small_components(i - 2,j)) + &
                 w2*(radial_small_components(i + 2,j) - radial_small_components(i - 1,j)) + &
                 w1*(radial_small_components(i + 1,j) - radial_small_components(i,j))
!
      end do
!
!   Set remaining deferred corrections to zero: slopes are
!   small in this region
!
      mfjm2 = orbital_last_grid_index(j) - 2
      dp(mfjm2:radial_point_count) = 0.0_real64
      dq(mfjm2:radial_point_count) = 0.0_real64
!
      return
      end subroutine defcor
