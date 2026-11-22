module read_specific_mass_shift_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_specific_mass_shift
contains

    ! WRITTEN  by C. Naz\'e  Oct. 2011
    !
    ! Translated by Gediminas Gaigalas 11/18/19
    subroutine read_specific_mass_shift(vinti_integrals,recoil_corrected_vinti_integrals)
        ! Module dependencies
        use calculation_file_units, only: two_body_cache_unit
        use grasp_parameters,    only: orbital_index_radix, max_orbitals
        use selected_level_indices, only: active_level_count
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use orbital_csf_basis, only: csf_count, max_orbitals
        use eigenstate_expansions, only: csf_mixing_coefficients
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: vinti_integrals(max_orbitals, max_orbitals), &
            & recoil_corrected_vinti_integrals(max_orbitals, max_orbitals)

        ! Local constants

        ! Local declarations
        real(real64) :: total_sms_contribution, leading_sms_contribution, coeffsms
        integer :: column_csf, row_csf, level_index, orbital_a, orbital_b, orbital_c, orbital_d, level_mixing_offset, &
            & packed_orbital_label, io_status

        call manage_angular_coefficient_buffer (1)

        rewind(two_body_cache_unit)
        do
            read (two_body_cache_unit,iostat = io_status) column_csf,row_csf
            if (io_status /= 0) exit
            ! cc        DO J = 1,NVEC
            read(two_body_cache_unit) coeffsms,packed_orbital_label
            orbital_d = mod (packed_orbital_label, orbital_index_radix)
            packed_orbital_label = packed_orbital_label/orbital_index_radix
            orbital_b  = mod (packed_orbital_label, orbital_index_radix)
            packed_orbital_label = packed_orbital_label/orbital_index_radix
            orbital_c = mod (packed_orbital_label, orbital_index_radix)
            orbital_a = packed_orbital_label/orbital_index_radix
            do level_index = 1,active_level_count
                level_mixing_offset = (level_index-1)*csf_count
                leading_sms_contribution = - csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)                     &
                    * coeffsms                                              &
                    * vinti_integrals (orbital_a,orbital_c)*vinti_integrals(orbital_b,orbital_d)
                total_sms_contribution = - csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)                       &
                    * coeffsms                                              &
                    * ( recoil_corrected_vinti_integrals(orbital_a,orbital_c)*vinti_integrals(orbital_b, &
                        & orbital_d)                        &
                    + recoil_corrected_vinti_integrals(orbital_b,orbital_d)*vinti_integrals(orbital_a,orbital_c))/2.0_real64
                if (row_csf/=column_csf) then
                    total_sms_contribution = 2.0_real64 * total_sms_contribution
                    leading_sms_contribution = 2.0_real64 * leading_sms_contribution
                end if
                leading_sms_by_level(level_index) = leading_sms_by_level(level_index) + leading_sms_contribution
                total_sms_by_level(level_index) = total_sms_by_level(level_index) + total_sms_contribution
            end do

        end do
        call manage_angular_coefficient_buffer (3)
        return
    end subroutine read_specific_mass_shift
end module read_specific_mass_shift_interface
