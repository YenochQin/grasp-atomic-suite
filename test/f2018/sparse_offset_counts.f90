program sparse_offset_counts
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use eigenstate_expansions, only: hamiltonian_energy_shift
    use rci_restart_records, only: scan_rci_matrix_block
    use sparse_hamiltonian_storage, only: stored_matrix_element_count
    implicit none (type, external)

    integer(int64), parameter :: max_default_count = int(huge(0),int64)
    integer :: unit, row_count, status, block_rows, row_index
    real(real64) :: last_shift

    ! Seed the running total and scan complete, valid records through production I/O.
    ! Zero, one and two rows leave totals immediately below, at and above the limit.
    do row_count = 0, 2
        open(newunit=unit, status="scratch", form="unformatted")
        write(unit) 2, 2, 0, 1
        do row_index = 1, row_count
            write(unit) 1, 0.0_real64, 1.0_real64, row_index
        end do
        rewind(unit)
        last_shift = 0.0_real64
        hamiltonian_energy_shift = 0.0_real64
        stored_matrix_element_count = max_default_count - 1_int64
        call scan_rci_matrix_block(unit, 1, 2, 2, 0, 1, status, block_rows, last_shift)
        if (status /= 0 .or. block_rows /= row_count) error stop "RCI count boundary fixture records"
        if (stored_matrix_element_count /= max_default_count - 1_int64 + row_count) &
            error stop "RCI accumulated count must retain all bits"
        close(unit)
    end do
end program sparse_offset_counts
