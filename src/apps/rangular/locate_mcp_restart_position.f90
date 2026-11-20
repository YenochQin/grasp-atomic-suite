!***********************************************************************
!                                                                      *
subroutine locate_mcp_restart_position(first_bra_csf_index, first_ket_csf_index, matrix_element_position, &
    & one_electron_coefficient_count, &
    & two_electron_counts_by_rank)
    ! Module dependencies
    use mcp_generation_options, only: mcp_temporary_units
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
    use orbital_csf_basis, only: csf_count
    implicit none (type, external)

    ! Arguments
    integer, intent(inout) :: first_bra_csf_index
    integer, intent(inout) :: first_ket_csf_index
    integer, intent(inout) :: matrix_element_position
    integer, intent(out) :: one_electron_coefficient_count
    integer, pointer, intent(in) :: two_electron_counts_by_rank(:)

    ! Local declarations
    integer :: k, io_status, icread, irread, indx, nrec, i, label
    integer, save :: sms_restart_unit = -1
    real(real64) :: coeff
    character :: srtlab*8, mcp_channel_label*3
    !-----------------------------------------------
    !
    !   Begin by examining file 30; this is the last to be updated by
    !   SUBROUTINE MCP
    !
    !   Read and check the character part of the file header
    !
    rewind (mcp_temporary_units(30))
    read (mcp_temporary_units(30)) mcp_channel_label, srtlab
    if (srtlab == "  SORTED") then
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
        read (mcp_temporary_units(30), iostat=io_status) icread, irread, indx
        if (io_status /= 0) exit
        first_bra_csf_index = icread
        first_ket_csf_index = irread
        matrix_element_position = indx
    end do

    ! All coefficients are present. Tell the caller to proceed to sorting.
    if (first_bra_csf_index == csf_count .and. first_ket_csf_index == csf_count) then
        first_bra_csf_index = csf_count + 1
        first_ket_csf_index = csf_count + 1
        return
    end if

    ! Position each coefficient channel after the records that will be kept.
    do k = 31, 32 + max_coulomb_multipole_rank
        rewind (mcp_temporary_units(k))
        nrec = 3
        do i = 1, nrec
            read (mcp_temporary_units(k))
        end do
        do
            read (mcp_temporary_units(k), iostat=io_status) indx, label, coeff
            if (io_status /= 0) exit
            if (indx > matrix_element_position) exit
            nrec = nrec + 1
        end do
        rewind (mcp_temporary_units(k))
        do i = 1, nrec
            read (mcp_temporary_units(k))
        end do
        if (k > 31) then
            two_electron_counts_by_rank(k - 32) = nrec - 3
        else
            one_electron_coefficient_count = nrec - 3
        end if
    end do

    ! The SMS file has the same number of data records as channel 33.
    if (sms_restart_unit == -1 .and. two_electron_counts_by_rank(1) > 0) then
        open(newunit=sms_restart_unit, file="fort.20", status="unknown", form="unformatted")
    end if
    if (sms_restart_unit /= -1) rewind (sms_restart_unit)
    do i = 1, two_electron_counts_by_rank(1)
        read (sms_restart_unit)
    end do

    first_ket_csf_index = first_ket_csf_index + 1
    if (first_ket_csf_index > csf_count) then
        first_bra_csf_index = first_bra_csf_index + 1
        first_ket_csf_index = first_bra_csf_index
    end if
end subroutine locate_mcp_restart_position
