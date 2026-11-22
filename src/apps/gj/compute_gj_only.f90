! Evaluate the diagonal atomic-state matrix elements for the Landé factor and its correction.
module compute_gj_only_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: compute_gj_only
contains

    subroutine compute_gj_only(lande_g_factors, lande_g_factor_corrections)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_orbitals
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation
        use eigenstate_expansions, only: csf_mixing_coefficients
        use first_order_csf_cutoffs,        only: first_order_csf_cutoff
        use standard_io_units,        only: stderr_unit
        use orbital_csf_basis,           only: csf_count, orbital_count
        use selected_level_indices,          only: active_level_count
        use integer_format, only: format_integer
        use csf_quantum_queries, only: csf_total_j_code

        ! Procedure interfaces
        use hyperfine_one_electron_element_interface, only: hyperfine_one_electron_element

        ! Module dependencies
        use one_electron_tensor_coefficients, only: calculate_one_electron_tensor_coefficients
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap

        ! Procedure interfaces
        use hyperfine_radial_integral_interface, only: hyperfine_radial_integral
        implicit none (type, external)

        ! Arguments
        real(real64), intent(out) :: lande_g_factors(active_level_count), lande_g_factor_corrections(active_level_count)

        ! Local constants
        real(real64), parameter :: contribution_magnitude_cutoff = 1.0e-10_real64

        ! Local declarations
        integer :: first_orbital_index, second_orbital_index, operator_parity, column_csf_index, text_length, row_csf_index, &
                & column_csf_j_code, row_csf_j_code, &
            & csf_j_code_difference, bra_orbital_index, &
            & ket_orbital_index, &
            level_index, mixing_offset
        real(real64) :: shell_coefficients(max_orbitals)
        real(real64) :: lande_g_radial_integrals(max_orbitals, max_orbitals), radial_overlap_integrals(max_orbitals, &
            & max_orbitals), &
            lande_g_angular_factors(max_orbitals, max_orbitals), lande_g_correction_angular_factors(max_orbitals, max_orbitals)
        real(real64) :: unused_hyperfine_angular_factor, lande_g_angular_factor, lande_g_correction_angular_factor, &
            & lande_g_csf_matrix_element, lande_g_correction_csf_matrix_element, lande_g_contribution, &
            lande_g_correction_contribution
        character :: number_text*11

        lande_g_factors(:active_level_count) = 0.0_real64
        lande_g_factor_corrections(:active_level_count) = 0.0_real64
        ! Precompute the radial integrals and angular factors required by
        ! the two one-body operators contributing to the Landé factor.
        do first_orbital_index = 1, orbital_count
            do second_orbital_index = 1, orbital_count
                lande_g_radial_integrals(first_orbital_index, &
                & second_orbital_index) = hyperfine_radial_integral(first_orbital_index,second_orbital_index,1)
                radial_overlap_integrals(first_orbital_index,second_orbital_index) = integrate_weighted_radial_overlap( &
                & first_orbital_index,second_orbital_index,0)
                call hyperfine_one_electron_element(first_orbital_index, 1, second_orbital_index, unused_hyperfine_angular_factor, &
                    & lande_g_angular_factor, lande_g_correction_angular_factor)
                lande_g_angular_factors(first_orbital_index,second_orbital_index) = lande_g_angular_factor
                lande_g_correction_angular_factors(first_orbital_index,second_orbital_index) = lande_g_correction_angular_factor
            end do
        end do
        operator_parity = 1
        ! Sum same-J CSF pairs that contribute to diagonal matrix elements in the atomic-state basis.
        do column_csf_index = 1, csf_count
            if (mod(column_csf_index,100) == 0) then
                call format_integer(column_csf_index, number_text, text_length)
                write (stderr_unit, *) "Column "//number_text(1:text_length)//" complete;"
            end if
            do row_csf_index = 1, csf_count
                if (use_first_order_csf_perturbation .and. column_csf_index>first_order_csf_cutoff .and. column_csf_index/= &
                & row_csf_index) cycle
                column_csf_j_code = csf_total_j_code(column_csf_index)
                row_csf_j_code = csf_total_j_code(row_csf_index)
                csf_j_code_difference = column_csf_j_code - row_csf_j_code
                if (.not.(csf_j_code_difference == 0 .and. row_csf_index >= column_csf_index)) cycle
                lande_g_csf_matrix_element = 0.0_real64
                lande_g_correction_csf_matrix_element = 0.0_real64
                call calculate_one_electron_tensor_coefficients(1, operator_parity, column_csf_index, row_csf_index, &
                & bra_orbital_index, &
                    & ket_orbital_index, shell_coefficients)
                if (bra_orbital_index == 0) cycle
                if (bra_orbital_index == ket_orbital_index) then
                    do first_orbital_index = 1, orbital_count
                        if (abs(shell_coefficients(first_orbital_index)) <= contribution_magnitude_cutoff) cycle
                        lande_g_csf_matrix_element = lande_g_csf_matrix_element + &
                            & lande_g_angular_factors(first_orbital_index,first_orbital_index) &
                            & *lande_g_radial_integrals(first_orbital_index, &
                            & first_orbital_index)*shell_coefficients(first_orbital_index)
                        lande_g_correction_csf_matrix_element = lande_g_correction_csf_matrix_element + &
                            & lande_g_correction_angular_factors(first_orbital_index,first_orbital_index) &
                            & *radial_overlap_integrals(first_orbital_index, &
                            & first_orbital_index)*shell_coefficients(first_orbital_index)
                    end do
                else
                    if (abs(shell_coefficients(1)) <= contribution_magnitude_cutoff) cycle
                    lande_g_csf_matrix_element = lande_g_csf_matrix_element + &
                        & lande_g_angular_factors(bra_orbital_index,ket_orbital_index) &
                            & *lande_g_radial_integrals(bra_orbital_index, &
                        & ket_orbital_index)*shell_coefficients(1)
                    lande_g_correction_csf_matrix_element = lande_g_correction_csf_matrix_element + &
                        & lande_g_correction_angular_factors(bra_orbital_index,ket_orbital_index) &
                            & *radial_overlap_integrals(bra_orbital_index, &
                        & ket_orbital_index)*shell_coefficients(1)
                end if
                if (abs(lande_g_csf_matrix_element) <= contribution_magnitude_cutoff .and. &
                    & abs(lande_g_correction_csf_matrix_element) <= contribution_magnitude_cutoff) cycle
                do level_index = 1, active_level_count
                    mixing_offset = (level_index - 1)*csf_count
                    if (row_csf_index /= column_csf_index) then
                        lande_g_contribution = lande_g_csf_matrix_element*(csf_mixing_coefficients(column_csf_index + &
                            & mixing_offset)*csf_mixing_coefficients(row_csf_index + mixing_offset) + &
                            & csf_mixing_coefficients(row_csf_index &
                            + mixing_offset)*csf_mixing_coefficients(column_csf_index + mixing_offset))
                        lande_g_correction_contribution = &
                            & lande_g_correction_csf_matrix_element*(csf_mixing_coefficients(column_csf_index + &
                            & mixing_offset)*csf_mixing_coefficients(row_csf_index + mixing_offset) + &
                            & csf_mixing_coefficients(row_csf_index&
                            + mixing_offset)*csf_mixing_coefficients(column_csf_index + mixing_offset))
                    else
                        lande_g_contribution = lande_g_csf_matrix_element*csf_mixing_coefficients(column_csf_index + &
                            & mixing_offset)*csf_mixing_coefficients(row_csf_index + mixing_offset)
                        lande_g_correction_contribution = &
                            & lande_g_correction_csf_matrix_element*csf_mixing_coefficients(column_csf_index + &
                            & mixing_offset)*csf_mixing_coefficients(row_csf_index + mixing_offset)
                    end if
                    lande_g_factors(level_index) = lande_g_factors(level_index) + lande_g_contribution
                    lande_g_factor_corrections(level_index) = lande_g_factor_corrections(level_index) + &
                & lande_g_correction_contribution
                end do
            end do
        end do
        return
    end subroutine compute_gj_only
end module compute_gj_only_interface
