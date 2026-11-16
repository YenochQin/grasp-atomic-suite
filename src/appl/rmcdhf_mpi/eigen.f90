!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function eigen (j)
!                                                                      *
!   This function computes an estimate of the energy of orbital J .    *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
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
      use orb_c
      use pote_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      use dpbdt_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: piece1, piece2, piece3, piece4, piece5
!-----------------------------------------------
!
!
!   Initialization
!
      mtp = orbital_last_grid_index(j)
!
!   Exchange term
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = (radial_large_components(i,j)*xq(i)-radial_small_components(i,j)*xp(i))*radial_derivative_over_radius(i)
      end do
      call quad (piece1)
      piece1 = speed_of_light_au*piece1
!
!   Direct term
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = (radial_large_components(i,j)**2 + radial_small_components(i,j)**2)*yp(i)*radial_derivative_over_radius(i)
      end do
      call quad (piece2)
!
!   Kinetic energy terms
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = radial_small_components(i,j)**2*radial_grid_derivative(i)
      end do
      call quad (piece3)
      piece3 = 2.0d00*speed_of_light_au*speed_of_light_au*piece3
!
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = (radial_large_components(i,j)*radial_small_components(i,j))*radial_derivative_over_radius(i)
      end do
      call quad (piece4)
      piece4 = -2.0d00*dble(orbital_kappa(j))*speed_of_light_au*piece4
!
      call dpbdt (j)
      ta(1) = 0.0d00
      do i = 2, mtp
         ta(i) = radial_large_components(i,j)*tb(i) - radial_small_components(i,j)*ta(i)
      end do
      call quad (piece5)
      piece5 = speed_of_light_au*piece5/transformed_grid_step
!
!   Assembly
!
      eigen = piece1 + piece2 + piece3 + piece4 + piece5
!
      return
!
      end function eigen
