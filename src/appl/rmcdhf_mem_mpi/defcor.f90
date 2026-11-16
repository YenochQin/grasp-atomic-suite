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
      use vast_kind_param, only:  double
      use def_c
      use grid_c
      use wave_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: w3 = 1.0d00/120.0d00
      real(double), parameter :: w2 = -15.0d00*w3
      real(double), parameter :: w1 = 40.0d00*w3
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
         dp(:2) = 0.0d00
         dq(:2) = 0.0d00
         first = .false.
      endif
!
!   Intermediate points
!
      mfjm3 = mf(j) - 3
      do i = 3, mfjm3
!
         dp(i) = w3*(pf(i + 3,j) - pf(i - 2,j)) + &
                 w2*(pf(i + 2,j) - pf(i - 1,j)) + &
                 w1*(pf(i + 1,j) - pf(i,j))
!
         dq(i) = w3*(qf(i + 3,j) - qf(i - 2,j)) + &
                 w2*(qf(i + 2,j) - qf(i - 1,j)) + &
                 w1*(qf(i + 1,j) - qf(i,j))
!
      end do
!
!   Set remaining deferred corrections to zero: slopes are
!   small in this region
!
      mfjm2 = mf(j) - 2
      dp(mfjm2:n) = 0.0d00
      dq(mfjm2:n) = 0.0d00
!
      return
      end subroutine defcor
