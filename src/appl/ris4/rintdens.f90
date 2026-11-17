!***********************************************************************
!                                                                      *
      real(real64) function rintdens (i, j)
!                                                                      *
!   The value of RINTDENS is an approximation to:                      *
!                                                                      *
!                                                                      *
!      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0              *
!                     I     J       I     J                            *
!                                                                      *
!   Call(s) to: [SMS92]:  POLINT                                       *
!                                                                      *
!   Written by Per Jonsson             Last revision: 24 Dec 1992      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c,           only: vacuum_speed_of_light_au, pi
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use wave_c, only: double, max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use polint_i, only: polint
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i, j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, ii, nmin
      real(real64) :: xa(3), ya(3)
      real(real64) :: dens
!-----------------------------------------------
!
      nmin = 10

      ii = 1
      do l = nmin,nmin+4,2
         xa(ii) = radial_grid(l)
         ya(ii) = (radial_large_components(l,i)*radial_large_components(l,j) + radial_small_components(l,i)* &
            & radial_small_components(l,j))/(4.0_real64*pi*radial_grid(l)*radial_grid(l))
         ii = ii + 1
      end do
      call polint (xa,ya,dens)
      rintdens = dens
!
      return
      end function rintdens
