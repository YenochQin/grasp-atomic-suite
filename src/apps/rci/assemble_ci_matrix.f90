! Assemble a CI matrix and record the sparse entries used for diagonalization.
module assemble_ci_matrix_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: assemble_ci_matrix
contains

    !        Generate Hamiltonian matrix for all blocks
    !        This routine calls setham to do the computation. It makes
    !        sure that the hamiltonian matrix is complete.
    !  Only node-0 has the correct, non-zero elsto. And in any case, elsto
    !  is not added to the matrix elements in files *.res
    !  This routine has been re-written to work in both single- and
    !  multi- processors.
    !
    ! Xinghong He 1998-06-23
    subroutine assemble_ci_matrix(inverse_nuclear_mass, symmetry_block_index, mpi_rank, mpi_process_count, &
        & closed_shell_breit_shift, restart_header_status, slf_en)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use eigenstate_expansions, only: hamiltonian_energy_shift
        use csf_block_cutoffs, only: active_csf_cutoffs
        use orbital_csf_basis, only: csf_count, orbital_count
        use ci_restart_units, only: restart_matrix_unit
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, block_stored_element_count, block_assembly_csf_count
        use standard_io_units, only: stderr_unit
        use file_positioning, only: position_calculation_file
        use rci_restart_records, only: scan_rci_matrix_block

        ! Procedure interfaces
        use assemble_hamiltonian_block_interface, only: assemble_hamiltonian_block
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: symmetry_block_index
        integer, intent(in) :: mpi_rank
        integer, intent(in) :: mpi_process_count
        integer, intent(out) :: restart_header_status
        real(real64)  :: inverse_nuclear_mass
        real(real64), intent(out) :: closed_shell_breit_shift
        real(real64) :: slf_en(*)

        ! Local declarations
        integer :: restart_rows_read, io_status, ipos, row_count, i, j, first_unrestored_csf_index

        stored_matrix_element_count = 0                      ! Counting continues in setham
        hamiltonian_energy_shift    = 0.0_real64
        closed_shell_breit_shift  = 0.0_real64

        call scan_rci_matrix_block(restart_matrix_unit, symmetry_block_index, csf_count, active_csf_cutoffs(1), &
            mpi_rank, mpi_process_count, io_status, restart_rows_read, closed_shell_breit_shift)
        ipos = 7 + 2*orbital_count + restart_rows_read
        if (io_status == 0) ipos = ipos + 1

        ! Find the maximum number of rows

        row_count = (csf_count - mpi_rank - 1 + mpi_process_count)/mpi_process_count
        if (csf_count < mpi_process_count) row_count = csf_count/(mpi_rank + 1)

        ! Report the number of rows read.
        ! A more suitable report on all nodes can be done here, but this will
        ! set a synchronization point.

        write (stderr_unit, *) restart_rows_read, " (total ", row_count, ") rows read from .res"
        if (mpi_rank == 0) write (summary_output_unit, *) restart_rows_read, " (total ", row_count, &
            ") rows read from .res"

        ! Position the file for the next record from setham

        do i = 1, symmetry_block_index - 1
            j = (csf_counts_by_block(i) - mpi_rank - 1 + mpi_process_count)/mpi_process_count
            if (csf_counts_by_block(i) < mpi_process_count) j = csf_counts_by_block(i)/(mpi_rank + 1)
            ipos = ipos + j + 1
        end do
        call position_calculation_file (0, restart_matrix_unit, ipos)

        if (io_status /= 0) write (restart_matrix_unit) csf_count, active_csf_cutoffs(1), mpi_rank, mpi_process_count

        if (restart_rows_read < row_count) then
            first_unrestored_csf_index = restart_rows_read*mpi_process_count + mpi_rank + 1
            ! ...Generate the rest of the Hamiltonian matrix
            call assemble_hamiltonian_block (mpi_rank, mpi_process_count, symmetry_block_index, closed_shell_breit_shift, &
                & first_unrestored_csf_index, &
                & stored_matrix_element_count, inverse_nuclear_mass, &
                slf_en)
        else
            block_stored_element_count = stored_matrix_element_count
            block_assembly_csf_count = csf_count
        end if

        restart_header_status = io_status

        return
    end subroutine assemble_ci_matrix
end module assemble_ci_matrix_interface
