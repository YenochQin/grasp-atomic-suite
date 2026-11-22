! Find a stored position in the MCP calculation files for restart handling.
module locate_mcp_restart_position_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: locate_mcp_restart_position
contains

    subroutine locate_mcp_restart_position(first_bra_csf_index, first_ket_csf_index, matrix_element_position, &
        & one_electron_coefficient_count, &
        & two_electron_counts_by_rank)
        ! Module dependencies
        use mcp_generation_options, only: mcp_temporary_units
        use, intrinsic :: iso_fortran_env, only: real64
        use mcp_generation_options, only: max_coulomb_multipole_rank
        use orbital_csf_basis, only: csf_count
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: first_bra_csf_index
        integer, intent(inout) :: first_ket_csf_index
        integer, intent(inout) :: matrix_element_position
        integer, intent(out) :: one_electron_coefficient_count
        integer, pointer :: two_electron_counts_by_rank(:)

        ! Local declarations
        integer :: channel_unit_index, io_status, stored_bra_csf_index, stored_ket_csf_index, &
            & stored_matrix_element_position, retained_record_count, record_index, packed_integral_label
        integer, save :: sms_restart_unit = -1
        real(real64) :: stored_coefficient
        character :: sort_status_label*8, mcp_channel_label*3
        ! Begin by examining file 30; this is the last to be updated by
        ! SUBROUTINE MCP
        !
        ! Read and check the character part of the file header
        rewind (mcp_temporary_units(30))
        read (mcp_temporary_units(30)) mcp_channel_label, sort_status_label
        if (sort_status_label == "  SORTED") then
            first_bra_csf_index = csf_count + 1
            first_ket_csf_index = csf_count + 1
            one_electron_coefficient_count = 0
            two_electron_counts_by_rank(:max_coulomb_multipole_rank) = 0
            return
        end if

        read (mcp_temporary_units(30))
        read (mcp_temporary_units(30))

        ! Retain the CSF pair and matrix position from the last complete record.
        do
            read (mcp_temporary_units(30), iostat=io_status) stored_bra_csf_index, stored_ket_csf_index, &
                & stored_matrix_element_position
            if (io_status /= 0) exit
            first_bra_csf_index = stored_bra_csf_index
            first_ket_csf_index = stored_ket_csf_index
            matrix_element_position = stored_matrix_element_position
        end do

        ! All coefficients are present. Tell the caller to proceed to sorting.
        if (first_bra_csf_index == csf_count .and. first_ket_csf_index == csf_count) then
            first_bra_csf_index = csf_count + 1
            first_ket_csf_index = csf_count + 1
            return
        end if

        ! Position each coefficient channel after the records that will be kept.
        do channel_unit_index = 31, 32 + max_coulomb_multipole_rank
            rewind (mcp_temporary_units(channel_unit_index))
            retained_record_count = 3
            do record_index = 1, retained_record_count
                read (mcp_temporary_units(channel_unit_index))
            end do
            do
                read (mcp_temporary_units(channel_unit_index), iostat=io_status) stored_matrix_element_position, &
                    & packed_integral_label, stored_coefficient
                if (io_status /= 0) exit
                if (stored_matrix_element_position > matrix_element_position) exit
                retained_record_count = retained_record_count + 1
            end do
            rewind (mcp_temporary_units(channel_unit_index))
            do record_index = 1, retained_record_count
                read (mcp_temporary_units(channel_unit_index))
            end do
            if (channel_unit_index > 31) then
                two_electron_counts_by_rank(channel_unit_index - 32) = retained_record_count - 3
            else
                one_electron_coefficient_count = retained_record_count - 3
            end if
        end do

        ! The SMS file has the same number of data records as channel 33.
        if (sms_restart_unit == -1 .and. two_electron_counts_by_rank(1) > 0) then
            open(newunit=sms_restart_unit, file="fort.20", status="unknown", form="unformatted")
        end if
        if (sms_restart_unit /= -1) rewind (sms_restart_unit)
        do record_index = 1, two_electron_counts_by_rank(1)
            read (sms_restart_unit)
        end do

        first_ket_csf_index = first_ket_csf_index + 1
        if (first_ket_csf_index > csf_count) then
            first_bra_csf_index = first_bra_csf_index + 1
            first_ket_csf_index = first_bra_csf_index
        end if
    end subroutine locate_mcp_restart_position
end module locate_mcp_restart_position_interface
