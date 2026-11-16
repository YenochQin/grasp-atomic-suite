!***********************************************************************
!                                                                      *
      subroutine setxuv(j)
!                                                                      *
!   This  SUBROUTINE  sets  up the arrays XU and XV, for use by  the   *
!   subprograms IN and  OUT.                                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
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
      use int_c
      use pote_c
      use scf_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, i
      real(real64) :: dmhh
!-----------------------------------------------
!
!   Define constants
!
      dmhh = -transformed_grid_step*0.5_real64
!
!   Set up arrays XU and XV; since XU(1), XV(1) are never used,
!   set them to some arbitrary value
!
      nm1 = radial_point_count - 1
      xu(1) = 0.0_real64
      xv(1) = 0.0_real64
      xu(2:nm1) = dmhh*(xp(3:nm1+1)*radial_derivative_over_radius(3:nm1+1)+xp(2:nm1)*radial_derivative_over_radius(2:nm1)) + dp(2&
         :nm1)
      xv(2:nm1) = dmhh*(xq(3:nm1+1)*radial_derivative_over_radius(3:nm1+1)+xq(2:nm1)*radial_derivative_over_radius(2:nm1)) + dq(2&
         :nm1)
!
      xu(radial_point_count) = dmhh*xp(radial_point_count)*radial_derivative_over_radius(radial_point_count) + &
         & dp(radial_point_count)
      xv(radial_point_count) = dmhh*xp(radial_point_count)*radial_derivative_over_radius(radial_point_count) + &
         & dq(radial_point_count)
!
      return
      end subroutine setxuv
