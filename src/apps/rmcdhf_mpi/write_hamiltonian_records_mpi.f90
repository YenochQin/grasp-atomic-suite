!***********************************************************************
!                                                                      *
subroutine write_hamiltonian_records_mpi(mpi_rank, mpi_process_count, block_csf_count)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   Routine for printing the Hamiltonian matrix.                       *
    !                                                                      *
    !   Written by Farid A Parpia             Last revision: 21 Dec 1992   *
    !   Block Version by Xinghong He          Last revision: 30 Jan 1999   *
    !                                                                      *
    !***********************************************************************
    use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
        matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
        nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
    !GG      USE mpi_C
    implicit none (type, external)

    ! Arguments
    integer  :: mpi_rank
    integer  :: mpi_process_count
    integer, intent(in) :: block_csf_count

    ! Local declarations
    integer :: ibeg, ico, idiag, list, iro
    !-----------------------------------------------
    !
    ibeg = 1
    do ico = mpi_rank + 1, block_csf_count, mpi_process_count
        idiag = stored_segment_end_offsets(ico)
        do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            call ensure_debug_output()
            write (angular_debug_unit,*) "H(",iro,ico,")= ", hamiltonian_values(list)
        end do
        ibeg = idiag + 1
    end do
    !
    return
end subroutine write_hamiltonian_records_mpi
