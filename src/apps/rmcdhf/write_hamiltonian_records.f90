!***********************************************************************
!                                                                      *
subroutine write_hamiltonian_records(mpi_rank, mpi_process_count, block_csf_count)
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
    implicit none (type, external)

    ! Arguments and local declarations
    integer, save :: first_column_debug_unit = -1
    integer, save :: second_column_debug_unit = -1
    integer, save :: third_column_debug_unit = -1
    integer  :: mpi_rank
    integer  :: mpi_process_count
    integer, intent(in) :: block_csf_count
    integer :: ibeg, ico, idiag, list, iro
    !-----------------------------------------------
    !
    !
    !  This is original, correct, full matrix
    !      ibeg = 1
    !      do ico = myid + 1, ncf, nprocs
    !         idiag = iendc(ico)
    !         do list = ibeg, idiag
    !            iro = irow(list)
    !            write (99,*) 'H(',iro,ico,')= ', emt(list)
    !         enddo
    !         ibeg = idiag + 1
    !      enddo
    !
    ibeg = 1
    do ico = 1, block_csf_count, 4
        ibeg = stored_segment_end_offsets(ico - 1) + 1
        idiag = stored_segment_end_offsets(ico)
        do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            if (first_column_debug_unit == -1) then
                open(newunit=first_column_debug_unit, file="fort.96", status="unknown", form="formatted")
            end if
            write (first_column_debug_unit, *) "H(", iro, ico, ")= ", hamiltonian_values(list)
        end do
        !ibeg = idiag + 1
    end do

    ibeg = stored_segment_end_offsets(1) + 1
    do ico = 2, block_csf_count, 4
        ibeg = stored_segment_end_offsets(ico - 1) + 1
        idiag = stored_segment_end_offsets(ico)
        do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            if (second_column_debug_unit == -1) then
                open(newunit=second_column_debug_unit, file="fort.97", status="unknown", form="formatted")
            end if
            write (second_column_debug_unit, *) "H(", iro, ico, ")= ", hamiltonian_values(list)
        end do
        !ibeg = idiag + 1
    end do

    ibeg = stored_segment_end_offsets(2) + 1
    do ico = 3, block_csf_count, 4
        ibeg = stored_segment_end_offsets(ico - 1) + 1
        idiag = stored_segment_end_offsets(ico)
        do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            if (third_column_debug_unit == -1) then
                open(newunit=third_column_debug_unit, file="fort.98", status="unknown", form="formatted")
            end if
            write (third_column_debug_unit, *) "H(", iro, ico, ")= ", hamiltonian_values(list)
        end do
        !ibeg = idiag + 1
    end do

    ibeg = stored_segment_end_offsets(3) + 1
    do ico = 4, block_csf_count, 4
        ibeg = stored_segment_end_offsets(ico - 1) + 1
        idiag = stored_segment_end_offsets(ico)
        do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            call ensure_debug_output()
            write (angular_debug_unit, *) "H(", iro, ico, ")= ", hamiltonian_values(list)
        end do
        !ibeg = idiag + 1
    end do
    return

end subroutine write_hamiltonian_records
