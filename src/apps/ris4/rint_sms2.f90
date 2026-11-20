!***********************************************************************
!                                                                      *
      real(real64) function rint_sms2 (first_orbital, second_orbital)
!                                                                      *
!   The value of RINT_SMS2 is an approximation to:                     *
!                                                                      *
!                                                                      *
!   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
!   Range.                                                             *
!                                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                                                      *
!   Written by  G. Gaigalas                                            *
!           and E. Gaudamauskas                                        *
!                        Last             revision:  09 October 2009   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use atomic_calculation_state, only: nuclear_charge, speed_of_light_au
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use sigma_1_i, only: sigma_1
      use sigma_2_i, only: sigma_2
      use radial_quadrature, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer ,intent(in)  :: first_orbital, second_orbital
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: l
      real(real64) :: result, apart1, apart2
!-----------------------------------------------
!
!   Tabulate integrand as required for SUBROUTINE QUAD
!
      mtp = min (orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
!
!   Value at first tabulation point is arbitrary
!
      call sigma_1(2,first_orbital,second_orbital,apart1)
      call sigma_2(2,first_orbital,second_orbital,apart2)
      ta(1) = 0.0_real64
      do 1 l = 2,mtp
!HFS        TA(L) = (R(L)**K)*(PF(L,I)*QF(L,J)+QF(L,I)*PF(L,J))*RP(L)
        ta(l) =  radial_grid_derivative(l)*(apart1*radial_small_components(l,first_orbital)*radial_large_components(l, &
           & second_orbital) &
           & &
              -         apart2*radial_large_components(l,first_orbital)*radial_small_components(l,second_orbital))/radial_grid(l)
    1 continue
!
!   Perform integration
!
      call quad (result)
      rint_sms2 = -result*nuclear_charge/speed_of_light_au
!
      return
      end function rint_sms2
