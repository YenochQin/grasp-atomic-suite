module read_seltzer_field_shift_factors_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_seltzer_field_shift_factors_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_seltzer_field_shift_factors_mpi(contact_density_integrals,total_nms_integrals,radius_integrals, &
        & radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals,&
        leading_nms_integrals,orbital_density_profiles,level_density_profiles,nuclear_density_point_count)
        ! Module dependencies
        use calculation_file_units, only: one_body_cache_unit
        !                                                                      *
        !   MPI version of DENSREAD_SELTZ. Each rank scans the full angular    *
        !   file but only accumulates records whose IC belongs to its strided   *
        !   subset.                                                             *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: orbital_index_radix, max_orbitals
        use debug_output_settings,          only: matrix_element_screening_cutoff
        use selected_level_indices, only: active_level_count
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use orbital_csf_basis, only: csf_count, max_orbitals
        use eigenstate_expansions, only: csf_mixing_coefficients
        use mpi_state,            only: mpi_rank, mpi_process_count
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nuclear_density_point_count
        real(real64), intent(in) :: contact_density_integrals(max_orbitals, max_orbitals), &
            & total_nms_integrals(max_orbitals, max_orbitals),  &
            radius_integrals(max_orbitals, max_orbitals), radius_squared_integrals(max_orbitals, max_orbitals), &
            inverse_radius_integrals(max_orbitals, max_orbitals), inverse_radius_squared_integrals(max_orbitals, max_orbitals),  &
            leading_nms_integrals(max_orbitals, max_orbitals)
        real(real64), intent(out) :: level_density_profiles(active_level_count, nuclear_density_point_count)
        real(real64), intent(in) :: orbital_density_profiles(max_orbitals, max_orbitals, nuclear_density_point_count)

        ! Local constants

        ! Local declarations
        real(real64) :: cached_shell_coefficients(max_orbitals)
        real(real64) :: level_density_profile_contribution(nuclear_density_point_count), &
            & csf_pair_density_profile(nuclear_density_point_count)
        real(real64) :: contact_density_element, total_nms_element, radius_element, radius_squared_element, inverse_radius_element
        real(real64) :: inverse_radius_squared_element, leading_nms_element
        real(real64) :: contact_density_contribution, total_nms_contribution, radius_contribution, &
            & radius_squared_contribution, inverse_radius_contribution
        real(real64) :: inverse_radius_squared_contribution, leading_nms_contribution
        real(real64) :: weighted_mixing_product, pair_scale
        integer :: io_status, ia, ib, column_csf, row_csf, i, level_index, l, level_mixing_offset, packed_orbital_label, &
            & coefficient_count
        logical :: rank_owns_record

        level_density_profiles(:,:) = 0.0_real64
        rewind (one_body_cache_unit)
        do
            read (one_body_cache_unit,iostat = io_status) column_csf,row_csf,coefficient_count
            if (io_status /= 0) exit
            rank_owns_record = mod(column_csf-1, mpi_process_count) == mpi_rank
            contact_density_element = 0.0_real64
            total_nms_element = 0.0_real64
            radius_element = 0.0_real64
            radius_squared_element = 0.0_real64
            inverse_radius_element = 0.0_real64
            inverse_radius_squared_element = 0.0_real64
            leading_nms_element = 0.0_real64
            csf_pair_density_profile(:) = 0.0_real64

            do i= 1,coefficient_count
                read(one_body_cache_unit) cached_shell_coefficients(i),packed_orbital_label
                if (.not. rank_owns_record) cycle
                ia = packed_orbital_label/orbital_index_radix
                ib = mod(packed_orbital_label,orbital_index_radix)
                contact_density_element = contact_density_element + contact_density_integrals(ia,ib)*cached_shell_coefficients(i)
                do l = 2,nuclear_density_point_count
                    csf_pair_density_profile(l) = csf_pair_density_profile(l) + orbital_density_profiles(ia,ib,l)* &
                        & cached_shell_coefficients(i)
                end do
                total_nms_element = total_nms_element + total_nms_integrals(ia,ib)*cached_shell_coefficients(i)
                radius_element = radius_element + radius_integrals(ia,ib)*cached_shell_coefficients(i)
                radius_squared_element = radius_squared_element + radius_squared_integrals(ia,ib)*cached_shell_coefficients(i)
                inverse_radius_element = inverse_radius_element + inverse_radius_integrals(ia,ib)*cached_shell_coefficients(i)
                inverse_radius_squared_element = inverse_radius_squared_element + inverse_radius_squared_integrals(ia,ib)* &
                    & cached_shell_coefficients(i)
                leading_nms_element = leading_nms_element + leading_nms_integrals(ia,ib)*cached_shell_coefficients(i)
            end do

            if (rank_owns_record) then
                if (abs(contact_density_element) <= matrix_element_screening_cutoff .and. &
                    abs(total_nms_element) <= matrix_element_screening_cutoff .and. &
                    abs(radius_element) <= matrix_element_screening_cutoff .and. &
                    abs(radius_squared_element) <= matrix_element_screening_cutoff .and. &
                    abs(inverse_radius_element) <= matrix_element_screening_cutoff .and. &
                    abs(inverse_radius_squared_element) <= matrix_element_screening_cutoff .and. &
                    abs(leading_nms_element) <= matrix_element_screening_cutoff .and. &
                    all(abs(csf_pair_density_profile(:)) <= matrix_element_screening_cutoff)) cycle
                if (row_csf == column_csf) then
                    pair_scale = 1.0_real64
                else
                    pair_scale = 2.0_real64
                end if
                do level_index = 1,active_level_count
                    level_mixing_offset = (level_index-1)*csf_count
                    weighted_mixing_product = pair_scale * csf_mixing_coefficients(column_csf+level_mixing_offset)* &
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

        return
    end subroutine read_seltzer_field_shift_factors_mpi
end module read_seltzer_field_shift_factors_mpi_interface
