!***********************************************************************
!                                                                      *
      subroutine maxarr(j)
!                                                                      *
!   This subroutine finds the least self-consistent orbital            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
!                                                                      *
! J initialized to zero
! XHH 1997.02.14
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use fixd_c, only: double, fixed_orbital_count, is_orbital_fixed, max_orbitals
      use orb_c
      use scf_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: dlrgst
!-----------------------------------------------
!
!
      j = 0
      dlrgst = 0.0_real64
      do i = 1, orbital_count
         if (is_orbital_fixed(i)) cycle
         if (scnsty(i) <= dlrgst) cycle
         dlrgst = scnsty(i)
         j = i
      end do
!
      return
      end subroutine maxarr
