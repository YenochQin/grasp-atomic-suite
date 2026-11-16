!***********************************************************************
!                                                                      *
      subroutine consis(j)
!                                                                      *
!   This routine computes the weighted self-consistency of orbital J   *
!                                                                      *
!   Written by Farid A Parpia, at OXFORD    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grid_c
      use int_c
      use scf_c
      use wave_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtp, i
      real(real64) :: scmea, deltao
!-----------------------------------------------
!
!
      scmea = 0.0_real64
      mtp = min(mtp0,orbital_last_grid_index(j))
      do i = 1, mtp
         deltao = abs(p(i)-radial_large_components(i,j)) + abs(q(i)-radial_small_components(i,j))
         scmea = dmax1(deltao,scmea)
      end do
      scnsty(j) = scmea*sqrt(ucf(j))
!
      return
      end subroutine consis
