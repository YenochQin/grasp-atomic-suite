module build_exchange_radial_potential_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_exchange_radial_potential
contains

    !***********************************************************************
    !                                                                      *
    subroutine build_exchange_radial_potential(orbital_index)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   This subroutine tabulates the exchange terms (the first terms on   *
        !   the  right-hand  sides of eqs (14), I P Grant, B J McKenzie, P H   *
        !   Norrington, D F Mayers, and N C Pyper, Computer  Phys  Commun 21   *
        !   (1980) 211) for orbital J. The exchange terms are stored  in the   *
        !   common arrays XP and XQ.                                           *
        !                                                                      *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: orbital_index_radix
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state, only: speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use scf_radial_potentials, only: large_exchange_source, small_exchange_source
        use scf_iteration_state, only: exchange_coefficient_count, exchange_coefficients, exchange_coefficient_labels
        use radial_quadrature_workspace, only: radial_workspace_b
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use hartree_yz_functions, only: evaluate_hartree_yz_functions
        use radial_ascii_plot, only: plot_radial_orbital_ascii
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index

        ! Local constants

        ! Local declarations
        integer :: i, coefficient_index, packed_coefficient_label, k, first_integral_orbital, second_integral_orbital, &
            source_orbital, nb2, row_count, ii, ii1, ii2, mfi
        real(real64) :: coeff, weighted_radial_integral
        !-----------------------------------------------
        !
        !   Debug printout: header
        !
        if (rmcdhf_debug_flags(27) .or. rmcdhf_debug_flags(28)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,/," Exchange potential contributions (coefficients will ",&
                &" be divided by C) for ",1i2,1a2," orbital :"/,/)') orbital_principal_n(orbital_index), &
                & orbital_symmetry_labels(orbital_index)
        end if
        !
        !   Clear for accumulation of sums
        !
        large_exchange_source(:radial_point_count) = 0.0_real64
        small_exchange_source(:radial_point_count) = 0.0_real64
        !
        !   Add contributions from exchange terms
        !
        do coefficient_index = 1, exchange_coefficient_count

            ! Decode information in label
            packed_coefficient_label = exchange_coefficient_labels(coefficient_index)
            k = mod(packed_coefficient_label,orbital_index_radix)
            packed_coefficient_label = packed_coefficient_label/orbital_index_radix
            first_integral_orbital = mod(packed_coefficient_label,orbital_index_radix)
            packed_coefficient_label = packed_coefficient_label/orbital_index_radix
            second_integral_orbital = mod(packed_coefficient_label,orbital_index_radix)
            source_orbital = packed_coefficient_label/orbital_index_radix
            coeff = exchange_coefficients(coefficient_index)

            ! Debug printout: composition
            if (rmcdhf_debug_flags(27)) then
                call ensure_debug_output()
                write (angular_debug_unit, '(/,25x,"(",1i2,")"/,1x,1p,d21.14,"* Y    (",1i2,1a2,",",1i2,1a2,") ","* P (",&
                    &1i2,1a2,")")') k, coeff, orbital_principal_n(first_integral_orbital), &
                        & orbital_symmetry_labels(first_integral_orbital), &
                    & orbital_principal_n(second_integral_orbital)&
                    , orbital_symmetry_labels(second_integral_orbital), orbital_principal_n(source_orbital), &
                        & orbital_symmetry_labels(source_orbital)
            end if

            call evaluate_hartree_yz_functions (k, first_integral_orbital, second_integral_orbital)
            !
            !   Accumulate contributions
            !
            mfi = orbital_last_grid_index(source_orbital)
            coeff = coeff/speed_of_light_au
            !DO I = 1, MF(IORB)
            do i = 1, mfi
                !GG         DO I = 1, N
                weighted_radial_integral = coeff*radial_workspace_b(i)
                large_exchange_source(i) = large_exchange_source(i) + &
                    & weighted_radial_integral*radial_small_components(i,source_orbital)
                small_exchange_source(i) = small_exchange_source(i) - &
                    & weighted_radial_integral*radial_large_components(i,source_orbital)
            end do
        end do
        !
        !   Debug printout: potential functions
        !
        if (rmcdhf_debug_flags(28)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,/,31x,"(P)",19x,"(Q)",41x,"(P)",19x,"(Q)"/,&
                &2(" --------- r --------- ------ X  (r) -------"," ------ X  (r) -------"))')
            nb2 = radial_point_count/2
            if (2*nb2 == radial_point_count) then
                row_count = nb2
            else
                row_count = nb2 + 1
            end if
            do ii = 1, row_count
                ii1 = ii
                ii2 = ii1 + row_count
                if (ii2 <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), large_exchange_source(ii1), &
                        small_exchange_source(ii1), radial_grid(ii2), large_exchange_source( &
                        & ii2), small_exchange_source(&
                        ii2)
                else if (ii1 <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), large_exchange_source(ii1), &
                        small_exchange_source(ii1)
                end if
            end do
            call plot_radial_orbital_ascii (large_exchange_source, 1.0_real64, small_exchange_source, speed_of_light_au, &
                & radial_point_count)
        end if
        !
        return
        !
        return
        !
    end subroutine build_exchange_radial_potential
end module build_exchange_radial_potential_interface
