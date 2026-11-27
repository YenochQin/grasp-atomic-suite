! Distribute the specific mass-shift matrix-element calculation across MPI ranks.
module calculate_specific_mass_shift_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_specific_mass_shift_mpi
contains

    ! MPI version of SMSNEW. Each rank handles a strided subset of the
    ! outer IC loop and can optionally write a per-rank angular-data
    ! fragment that the root process later merges into name.isotope.two-body.angular.cache.
    subroutine calculate_specific_mass_shift_mpi(angular_cache_write_mode,vinti_integrals,recoil_corrected_vinti_integrals)
        ! Module dependencies
        use calculation_file_units, only: two_body_cache_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: orbital_index_radix, max_orbitals
        use debug_output_settings, only: matrix_element_screening_cutoff
        use selected_level_indices, only: active_level_count
        use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count
        use angular_coefficient_buffer_state, only: angular_coefficient_count, angular_coefficient_labels, angular_coefficients
        use eigenstate_expansions, only: csf_mixing_coefficients
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use mpi_state,            only: mpi_rank, mpi_process_count
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
        use csf_quantum_queries, only: csf_parity_code, csf_total_j_code
        use integer_format, only: format_integer
        use two_electron_angular_coefficients, only: calculate_two_electron_angular_coefficients
        use closed_shell_mcp, only: calculate_closed_shell_mcp
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: vinti_integrals(max_orbitals, max_orbitals), &
            & recoil_corrected_vinti_integrals(max_orbitals, max_orbitals)
        integer, intent(in) :: angular_cache_write_mode

        ! Local constants

        ! Local declarations
        character(len=11) :: number_text
        real(real64) :: weighted_mixing_product(active_level_count)
        real(real64) :: total_sms_contribution, leading_sms_contribution, angular_coefficient, radial_k1, radial_sms,  &
            pair_scale, term_k1, term_sms
        integer :: level_index,k,level_mixing_offset,csf_index,row_csf,column_two_j_plus_one,include_closed_shell_terms, &
            & text_length, orbital_a, orbital_b, orbital_c, orbital_d
        integer :: orbital_index, column_csf_occupation, row_csf_occupation, differing_subshell_count

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
                        if (iabs(column_csf_occupation - row_csf_occupation) > 2) then
                            differing_subshell_count = 5
                            exit
                        end if
                        differing_subshell_count = differing_subshell_count + iabs(column_csf_occupation - row_csf_occupation)
                        if (differing_subshell_count > 4) exit
                    end do
                    if (differing_subshell_count > 4 .or. mod(differing_subshell_count,2) /= 0) cycle

                    angular_coefficient_count = 0
                    call calculate_two_electron_angular_coefficients (csf_index, row_csf, calculate_closed_shell_mcp, &
                        & include_closed_shell_terms, 1)
                    if (angular_coefficient_count == 0) cycle
                    if (row_csf == csf_index) then
                        pair_scale = 1.0_real64
                    else
                        pair_scale = 2.0_real64
                    end if
                    do level_index = 1,active_level_count
                        level_mixing_offset = (level_index-1)*csf_count
                        weighted_mixing_product(level_index) = pair_scale * csf_mixing_coefficients(csf_index+ &
                            & level_mixing_offset) * csf_mixing_coefficients(row_csf+level_mixing_offset)
                    end do
                    do k = 1,angular_coefficient_count
                        angular_coefficient = angular_coefficients(k)
                        if (abs (angular_coefficient) > matrix_element_screening_cutoff) then
                            if (angular_coefficient_labels(5,k) == 1) then
                                orbital_a = angular_coefficient_labels(1,k)
                                orbital_b = angular_coefficient_labels(2,k)
                                orbital_c = angular_coefficient_labels(3,k)
                                orbital_d = angular_coefficient_labels(4,k)
                                radial_k1 = vinti_integrals(orbital_a,orbital_c) * vinti_integrals(orbital_b,orbital_d)
                                radial_sms = (recoil_corrected_vinti_integrals(orbital_a,orbital_c) * &
                                    & vinti_integrals(orbital_b,orbital_d) +   &
                                    recoil_corrected_vinti_integrals(orbital_b,orbital_d) * vinti_integrals(orbital_a, &
                                        & orbital_c))/2.0_real64
                                if (radial_k1 == 0.0_real64 .and.                     &
                                    radial_sms == 0.0_real64) cycle
                                term_k1 = -angular_coefficient * radial_k1
                                term_sms = -angular_coefficient * radial_sms
                                if(angular_cache_write_mode==1) write(two_body_cache_unit) csf_index,row_csf
                                if(angular_cache_write_mode==1) then
                                    write(two_body_cache_unit) angular_coefficient, ((orbital_a*orbital_index_radix + &
                                        & orbital_c)*orbital_index_radix+orbital_b)* &
                                        & orbital_index_radix+orbital_d
                                end if
                                do level_index = 1,active_level_count
                                    leading_sms_contribution = weighted_mixing_product(level_index) * term_k1
                                    total_sms_contribution = weighted_mixing_product(level_index) * term_sms
                                    leading_sms_by_level(level_index) = leading_sms_by_level(level_index) + leading_sms_contribution
                                    total_sms_by_level(level_index) = total_sms_by_level(level_index) + total_sms_contribution
                                end do
                            end if
                        end if
                    end do
                end if
            end do
        end do

        call manage_angular_coefficient_buffer (3)
        return
    end subroutine calculate_specific_mass_shift_mpi
end module calculate_specific_mass_shift_mpi_interface
