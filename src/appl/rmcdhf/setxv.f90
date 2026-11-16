!***********************************************************************
!                                                                      *
      subroutine setxv(j)
!                                                                      *
!   This  subprogram  sets up the inhomogeneous terms for the varia-   *
!   tion equations.                                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use grid_c
      use int_c
      use wave_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: hhc
!-----------------------------------------------
!
      hhc = 0.5_real64*transformed_grid_step/speed_of_light_au
!
!   Set up arrays TF and TG
!
      do i = 1, radial_point_count
         xu(i) = -radial_small_components(i,j)*hhc*radial_grid_derivative(i)
         xv(i) = radial_large_components(i,j)*hhc*radial_grid_derivative(i)
      end do
!
      xu(:radial_point_count-1) = xu(2:radial_point_count) + xu(:radial_point_count-1)
      xv(:radial_point_count-1) = xv(2:radial_point_count) + xv(:radial_point_count-1)
!
      return
!
      end subroutine setxv
