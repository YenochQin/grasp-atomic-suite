! Write the Hamiltonian records assigned to this MPI rank.
module write_hamiltonian_records_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_hamiltonian_records_mpi
contains

    ! Routine for printing the Hamiltonian matrix.
    !
    ! Written by Farid A Parpia. Last revision: 21 Dec 1992
    ! Block Version by Xinghong He. Last revision: 30 Jan 1999
    subroutine write_hamiltonian_records_mpi(mpi_rank, mpi_process_count, block_csf_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: int64
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use sparse_hamiltonian_storage, only: hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices
        implicit none (type, external)

        ! Arguments
        integer  :: mpi_rank
        integer  :: mpi_process_count
        integer, intent(in) :: block_csf_count

        ! Local declarations
        integer :: column_index, row_index
        integer(int64) :: first_matrix_element_index, diagonal_element_index, matrix_element_index

        first_matrix_element_index = 1
        do column_index = mpi_rank + 1, block_csf_count, mpi_process_count
            diagonal_element_index = stored_segment_end_offsets(column_index)
            do matrix_element_index = first_matrix_element_index, diagonal_element_index
                row_index = matrix_partner_indices(matrix_element_index)
                call ensure_debug_output()
                write (angular_debug_unit,*) "H(",row_index,column_index,")= ", hamiltonian_values(matrix_element_index)
            end do
            first_matrix_element_index = diagonal_element_index + 1
        end do
        return
    end subroutine write_hamiltonian_records_mpi
end module write_hamiltonian_records_mpi_interface
