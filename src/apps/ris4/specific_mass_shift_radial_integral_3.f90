module specific_mass_shift_radial_integral_3_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: specific_mass_shift_radial_integral_3
contains

    !***********************************************************************
    !                                                                      *
    real(real64) function specific_mass_shift_radial_integral_3 (first_orbital, second_orbital)
        ! Module dependencies
        !                                                                      *
        !   The value of RINT_SMS3 is an approximation to:                     *
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
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_grid, radial_grid_derivative
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use atomic_calculation_state, only: nuclear_charge, speed_of_light_au
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital, second_orbital

        ! Local declarations
        integer      :: l
        real(real64) :: result
        !-----------------------------------------------
        !
        !   Tabulate integrand as required for SUBROUTINE QUAD
        !
        integration_last_index = min (orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        !
        !   Value at first tabulation point is arbitrary
        !
        radial_workspace_a(1) = 0.0_real64
        do l = 2,integration_last_index
            !HFS        TA(L) = (R(L)**K)*(PF(L,I)*QF(L,J)+QF(L,I)*PF(L,J))*RP(L)
            radial_workspace_a(l) =  radial_grid_derivative(l)*(radial_small_components(l, &
                first_orbital)*radial_large_components(l,second_orbital) &
                -         radial_large_components(l,first_orbital)*radial_small_components(l,second_orbital))/radial_grid(l)
        end do
        !
        !   Perform integration
        !
        call integrate_radial_grid (result)
        specific_mass_shift_radial_integral_3 = result*nuclear_charge/speed_of_light_au
        !
        return
    end function specific_mass_shift_radial_integral_3
end module specific_mass_shift_radial_integral_3_interface
