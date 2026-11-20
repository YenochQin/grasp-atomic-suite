!*******************************************************************
!                                                                  *
subroutine read_mcp_coefficients_from_memory(jblock)
    ! Module dependencies
    !                                                                  *
    !   Written by  G. Gaigalas               Vilnius, September 2021  *
    !                                                                  *
    !*******************************************************************
    !
    use sparse_hamiltonian_storage,     only: matrix_partner_indices, stored_segment_end_offsets
    use rmcdhf_mcp_memory, only: iendc_30, imax1_30, imax2_30, imin1_30, imin2_30, irow_30
    use serial_mpi_context, only: mpi_process_count, mpi_rank
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: jblock

    ! Local declarations
    integer :: i
    !-----------------------------------------------

    do i = imin1_30(jblock)+1+mpi_rank, imax1_30(jblock), mpi_process_count
        stored_segment_end_offsets(i-imin1_30(jblock)) = iendc_30(i)
    end do
    do i = imin2_30(jblock)+1, imax2_30(jblock)
        matrix_partner_indices(i-imin2_30(jblock)) = irow_30(i)
    end do

    return
end subroutine read_mcp_coefficients_from_memory
