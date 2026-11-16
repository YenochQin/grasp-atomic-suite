!***********************************************************************
!                                                                      *
      subroutine setxz(j)
!                                                                      *
!   This subprogram sets the inhomogeneous terms to zero.              *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grid_c
      use int_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!
      xu(:radial_point_count) = 0.0_real64
      xv(:radial_point_count) = 0.0_real64
!
      return
      end subroutine setxz
