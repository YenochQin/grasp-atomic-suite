module read_specific_mass_shift_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_specific_mass_shift_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_specific_mass_shift_mpi(vinti_integrals,recoil_corrected_vinti_integrals)
        ! Module dependencies
        use calculation_file_units, only: two_body_cache_unit
        !                                                                      *
        !   MPI version of SMSREAD. Each rank scans the full angular file but  *
        !   only accumulates records whose IC belongs to its strided subset.    *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: orbital_index_radix, max_orbitals
        use selected_level_indices, only: active_level_count
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use orbital_csf_basis, only: csf_count, max_orbitals
        use eigenstate_expansions, only: csf_mixing_coefficients
        use mpi_state,            only: mpi_rank, mpi_process_count
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: vinti_integrals(max_orbitals, max_orbitals), &
            & recoil_corrected_vinti_integrals(max_orbitals, max_orbitals)

        ! Local constants

        ! Local declarations
        real(real64) :: weighted_mixing_product(active_level_count)
        real(real64) :: total_sms_contribution, leading_sms_contribution, coeffsms, radial_k1, radial_sms,&
            pair_scale, term_k1, term_sms
        integer :: column_csf, row_csf, level_index, orbital_a, orbital_b, orbital_c, orbital_d, level_mixing_offset, &
            & packed_orbital_label, io_status
        logical :: rank_owns_record

        call manage_angular_coefficient_buffer (1)
        rewind(two_body_cache_unit)
        do
            read (two_body_cache_unit,iostat = io_status) column_csf,row_csf
            if (io_status /= 0) exit
            read(two_body_cache_unit) coeffsms,packed_orbital_label
            rank_owns_record = mod(column_csf-1, mpi_process_count) == mpi_rank
            if (rank_owns_record) then
                orbital_d = mod (packed_orbital_label, orbital_index_radix)
                packed_orbital_label = packed_orbital_label/orbital_index_radix
                orbital_b  = mod (packed_orbital_label, orbital_index_radix)
                packed_orbital_label = packed_orbital_label/orbital_index_radix
                orbital_c = mod (packed_orbital_label, orbital_index_radix)
                orbital_a = packed_orbital_label/orbital_index_radix
                radial_k1 = vinti_integrals(orbital_a,orbital_c) * vinti_integrals(orbital_b,orbital_d)
                radial_sms = (recoil_corrected_vinti_integrals(orbital_a,orbital_c) * vinti_integrals(orbital_b, &
                    & orbital_d)                 &
                    + recoil_corrected_vinti_integrals(orbital_b,orbital_d) * vinti_integrals(orbital_a,orbital_c))/2.0_real64
                if (radial_k1 == 0.0_real64 .and. radial_sms == 0.0_real64) cycle
                term_k1 = -coeffsms * radial_k1
                term_sms = -coeffsms * radial_sms
                if (row_csf == column_csf) then
                    pair_scale = 1.0_real64
                else
                    pair_scale = 2.0_real64
                end if
                do level_index = 1,active_level_count
                    level_mixing_offset = (level_index-1)*csf_count
                    weighted_mixing_product(level_index) = pair_scale * csf_mixing_coefficients(column_csf+ &
                        & level_mixing_offset) * csf_mixing_coefficients(row_csf+level_mixing_offset)
                end do
                do level_index = 1,active_level_count
                    leading_sms_contribution = weighted_mixing_product(level_index) * term_k1
                    total_sms_contribution = weighted_mixing_product(level_index) * term_sms
                    leading_sms_by_level(level_index) = leading_sms_by_level(level_index) + leading_sms_contribution
                    total_sms_by_level(level_index) = total_sms_by_level(level_index) + total_sms_contribution
                end do
            end if
        end do
        call manage_angular_coefficient_buffer (3)
        return
    end subroutine read_specific_mass_shift_mpi
end module read_specific_mass_shift_mpi_interface
