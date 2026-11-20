!***********************************************************************
!                                                                      *
      subroutine rintdensvec (i,j,dint1vec,nrnuc)
!                                                                      *
!   The value of DINT1VEC is an approximation to:                      *
!                                                                      *
!                                                                      *
!      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) |                     *
!                       I     J       I     J                          *
!                                                                      *
!                                                                      *
!                                                                      *
!   Written by Per Jonsson  and Jörgen Ekman                           *
!   Last revision: 24 Dec 2013                                         *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_orbitals
      use atomic_calculation_state,            only: pi
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(out) :: dint1vec(max_orbitals, max_orbitals, radial_point_count)
      integer, intent(in) :: i, j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l
      integer, intent(in) :: nrnuc
!-----------------------------------------------
!
      do l = 2, nrnuc
         dint1vec(i,j,l) = (radial_large_components(l,i)*radial_large_components(l,j)+radial_small_components(l,i)* &
            & radial_small_components(l,j))/          &
                   (4.0_real64*pi*radial_grid(l)*radial_grid(l))
      end do
!
      return
      end subroutine rintdensvec
