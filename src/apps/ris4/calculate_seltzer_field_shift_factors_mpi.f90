module calculate_seltzer_field_shift_factors_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_seltzer_field_shift_factors_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_seltzer_field_shift_factors_mpi(angular_cache_write_mode,contact_density_integrals,total_nms_integrals, &
        & radius_integrals,radius_squared_integrals,inverse_radius_integrals, &
        inverse_radius_squared_integrals,leading_nms_integrals,orbital_density_profiles,level_density_profiles, &
            & nuclear_density_point_count)
        ! Module dependencies
        !                                                                      *
        !   MPI version of DENSNEW_SELTZ. Each rank handles a strided subset   *
        !   of the outer IC loop. In MPI mode angular coefficients are         *
        !   evaluated on the fly and not saved.                                *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: orbital_index_radix, max_orbitals
        use debug_output_settings, only: matrix_element_screening_cutoff
        use eigenstate_expansions, only: csf_mixing_coefficients
        use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count, orbital_kappa
        use selected_level_indices, only: active_level_count
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use mpi_state,            only: mpi_rank, mpi_process_count
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
        use integer_format, only: format_integer
        use csf_quantum_queries, only: csf_parity_code, csf_total_j_code
        use one_electron_scalar_coefficients, only: calculate_one_electron_scalar_coefficients
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: angular_cache_write_mode, nuclear_density_point_count
        real(real64), intent(in) :: contact_density_integrals(max_orbitals, max_orbitals), &
            & total_nms_integrals(max_orbitals, max_orbitals), &
            radius_integrals(max_orbitals, max_orbitals), radius_squared_integrals(max_orbitals, max_orbitals), &
            inverse_radius_integrals(max_orbitals, max_orbitals), inverse_radius_squared_integrals(max_orbitals, max_orbitals), &
            leading_nms_integrals(max_orbitals, max_orbitals)
        real(real64), intent(out) :: level_density_profiles(active_level_count, nuclear_density_point_count)
        real(real64), intent(in) :: orbital_density_profiles(max_orbitals, max_orbitals, nuclear_density_point_count)

        ! Local constants

        ! Local declarations
        real(real64) :: shell_coefficients(max_orbitals), stored_shell_coefficients(max_orbitals)
        real(real64) :: level_density_profile_contribution(nuclear_density_point_count), &
            & csf_pair_density_profile(nuclear_density_point_count)
        real(real64) :: contact_density_element, total_nms_element, radius_element, radius_squared_element, &
            & inverse_radius_element,     &
            inverse_radius_squared_element, leading_nms_element, contact_density_contribution, total_nms_contribution, &
                & radius_contribution,     &
            radius_squared_contribution, inverse_radius_contribution, inverse_radius_squared_contribution, &
                & leading_nms_contribution, weighted_mixing_product,      &
            pair_scale
        logical :: has_contr
        character :: number_text*11
        integer :: contributing_orbitals(max_orbitals)
        integer :: operator_rank, operator_parity, include_closed_shell_terms, csf_index, text_length, &
            & column_two_j_plus_one, row_csf, ia, &
            & ket_orbital_index, &
            & orbital_index, level_index
        integer :: l, level_mixing_offset, coefficient_count, column_csf_occupation, row_csf_occupation, differing_subshell_count

        level_density_profiles(:,:) = 0.0_real64
        operator_rank = 0
        operator_parity = 1
        include_closed_shell_terms = 1

        call manage_angular_coefficient_buffer (1)

        do csf_index = mpi_rank + 1, csf_count, mpi_process_count
            if (mpi_rank == 0 .and. mod(csf_index,100) == 1) then
                call format_integer (csf_index,number_text,text_length)
                print *, "Column "//number_text(1:text_length)//" complete on root;"
            end if

            column_two_j_plus_one = csf_total_j_code (csf_index)
            do row_csf = csf_index,csf_count
                if (csf_total_j_code(row_csf) == column_two_j_plus_one) then
                    if (csf_parity_code(csf_index) /= csf_parity_code(row_csf)) cycle
                    differing_subshell_count = 0
                    do orbital_index = 1,orbital_count
                        column_csf_occupation = csf_subshell_occupations(orbital_index,csf_index)
                        row_csf_occupation = csf_subshell_occupations(orbital_index,row_csf)
                        if (column_csf_occupation == row_csf_occupation) cycle
                        if (iabs(column_csf_occupation - row_csf_occupation) > 1) then
                            differing_subshell_count = 3
                            exit
                        end if
                        differing_subshell_count = differing_subshell_count + 1
                        if (differing_subshell_count > 2) exit
                    end do
                    if (differing_subshell_count /= 0 .and. differing_subshell_count /= 2) cycle
                    if (differing_subshell_count == 0 .and. csf_index /= row_csf) cycle

                    contact_density_element = 0.0_real64
                    total_nms_element = 0.0_real64
                    radius_element = 0.0_real64
                    radius_squared_element = 0.0_real64
                    inverse_radius_element = 0.0_real64
                    inverse_radius_squared_element = 0.0_real64
                    leading_nms_element = 0.0_real64
                    csf_pair_density_profile(:) = 0.0_real64
                    has_contr = .false.

                    call calculate_one_electron_scalar_coefficients(csf_index,row_csf,ia,ket_orbital_index,shell_coefficients)
                    if (ia /= 0) then
                        if (ia == ket_orbital_index) then
                            coefficient_count = 0
                            do ia = 1,orbital_count
                                if (abs (shell_coefficients(ia)) > matrix_element_screening_cutoff) then
                                    coefficient_count = coefficient_count + 1
                                    stored_shell_coefficients(coefficient_count) = shell_coefficients(ia)
                                    contributing_orbitals(coefficient_count) = ia
                                    contact_density_element = contact_density_element + contact_density_integrals(ia,ia)* &
                                        & shell_coefficients(ia)
                                    do l = 2,nuclear_density_point_count
                                        csf_pair_density_profile(l)=csf_pair_density_profile(l)+orbital_density_profiles(ia,ia, &
                                            & l)*shell_coefficients(ia)
                                    end do
                                    total_nms_element = total_nms_element + total_nms_integrals(ia,ia)*shell_coefficients(ia)
                                    radius_element = radius_element + radius_integrals(ia,ia)*shell_coefficients(ia)
                                    radius_squared_element = radius_squared_element + radius_squared_integrals(ia,ia)* &
                                        & shell_coefficients(ia)
                                    inverse_radius_element = inverse_radius_element + inverse_radius_integrals(ia,ia)* &
                                        & shell_coefficients(ia)
                                    inverse_radius_squared_element = inverse_radius_squared_element + &
                                        & inverse_radius_squared_integrals(ia,ia)*shell_coefficients(ia)
                                    leading_nms_element = leading_nms_element + leading_nms_integrals(ia,ia)*shell_coefficients(ia)
                                    has_contr = .true.
                                end if
                            end do
                        else
                            if (abs (shell_coefficients(1)) > matrix_element_screening_cutoff) then
                                if (orbital_kappa(ia)==orbital_kappa(ket_orbital_index)) then
                                    contact_density_element = contact_density_element + contact_density_integrals(ia, &
                                        & ket_orbital_index)*shell_coefficients(1)
                                    do l = 2,nuclear_density_point_count
                                        csf_pair_density_profile(l)=csf_pair_density_profile(l)+orbital_density_profiles(ia, &
                                            & ket_orbital_index,l)* &
                                            & shell_coefficients(1)
                                    end do
                                    total_nms_element = total_nms_element + total_nms_integrals(ia,ket_orbital_index)* &
                                        & shell_coefficients(1)
                                    radius_element = radius_element + radius_integrals(ia,ket_orbital_index)*shell_coefficients(1)
                                    radius_squared_element = radius_squared_element + radius_squared_integrals(ia, &
                                        & ket_orbital_index)*shell_coefficients(1)
                                    inverse_radius_element = inverse_radius_element + inverse_radius_integrals(ia, &
                                        & ket_orbital_index)*shell_coefficients(1)
                                    inverse_radius_squared_element = inverse_radius_squared_element + &
                                        & inverse_radius_squared_integrals(ia,ket_orbital_index)*shell_coefficients(1)
                                    leading_nms_element = leading_nms_element + leading_nms_integrals(ia, &
                                        & ket_orbital_index)*shell_coefficients(1)
                                    has_contr = .true.
                                end if
                            end if
                        end if
                    end if
                    if (.not. has_contr) cycle
                    if (abs(contact_density_element) <= matrix_element_screening_cutoff .and. &
                        abs(total_nms_element) <= matrix_element_screening_cutoff .and. &
                        abs(radius_element) <= matrix_element_screening_cutoff .and. &
                        abs(radius_squared_element) <= matrix_element_screening_cutoff .and. &
                        abs(inverse_radius_element) <= matrix_element_screening_cutoff .and. &
                        abs(inverse_radius_squared_element) <= matrix_element_screening_cutoff .and. &
                        abs(leading_nms_element) <= matrix_element_screening_cutoff .and. &
                        all(abs(csf_pair_density_profile(:)) <= matrix_element_screening_cutoff)) cycle
                    if (row_csf == csf_index) then
                        pair_scale = 1.0_real64
                    else
                        pair_scale = 2.0_real64
                    end if
                    do level_index = 1,active_level_count
                        level_mixing_offset = (level_index-1)*csf_count
                        weighted_mixing_product = pair_scale * csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+level_mixing_offset)
                        contact_density_contribution = weighted_mixing_product*contact_density_element
                        level_density_profile_contribution(:) = weighted_mixing_product*csf_pair_density_profile(:)
                        total_nms_contribution = weighted_mixing_product*total_nms_element
                        radius_contribution = weighted_mixing_product*radius_element
                        radius_squared_contribution = weighted_mixing_product*radius_squared_element
                        inverse_radius_contribution = weighted_mixing_product*inverse_radius_element
                        inverse_radius_squared_contribution = weighted_mixing_product*inverse_radius_squared_element
                        leading_nms_contribution = weighted_mixing_product*leading_nms_element
                        contact_density_by_level(level_index) = contact_density_by_level(level_index) + contact_density_contribution
                        do l = 2,nuclear_density_point_count
                            level_density_profiles(level_index,l) = level_density_profiles(level_index,l) + &
                                & level_density_profile_contribution(l)
                        end do
                        total_nms_by_level(level_index) = total_nms_by_level(level_index) + total_nms_contribution
                        radius_by_level(level_index) = radius_by_level(level_index) + radius_contribution
                        radius_squared_by_level(level_index) = radius_squared_by_level(level_index) + radius_squared_contribution
                        inverse_radius_by_level(level_index) = inverse_radius_by_level(level_index) + inverse_radius_contribution
                        inverse_radius_squared_by_level(level_index) = inverse_radius_squared_by_level(level_index) + &
                            & inverse_radius_squared_contribution
                        leading_nms_by_level(level_index) = leading_nms_by_level(level_index) + leading_nms_contribution
                    end do
                end if
            end do
        end do

        call manage_angular_coefficient_buffer (3)
        return
    end subroutine calculate_seltzer_field_shift_factors_mpi
end module calculate_seltzer_field_shift_factors_mpi_interface
