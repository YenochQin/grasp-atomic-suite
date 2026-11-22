! Write sparse Hamiltonian elements and indices for later diagonalization or restart.
module write_hamiltonian_records_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_hamiltonian_records
contains

    ! Routine for printing the Hamiltonian matrix.
    !
    ! Written by Farid A Parpia. Last revision: 21 Dec 1992
    ! Block Version by Xinghong He. Last revision: 30 Jan 1999
    subroutine write_hamiltonian_records(mpi_rank, mpi_process_count, block_csf_count)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use sparse_hamiltonian_storage, only: hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices
        implicit none (type, external)

        ! Arguments and local declarations
        integer, save :: first_column_debug_unit = -1
        integer, save :: second_column_debug_unit = -1
        integer, save :: third_column_debug_unit = -1
        integer  :: mpi_rank
        integer  :: mpi_process_count
        integer, intent(in) :: block_csf_count
        integer :: first_matrix_element_index, column_index, diagonal_element_index, matrix_element_index, row_index
        ! This is original, correct, full matrix
        !     ibeg = 1
        !        idiag = iendc(ico)
        !           iro = irow(list)
        !        ibeg = idiag + 1
        first_matrix_element_index = 1
        do column_index = 1, block_csf_count, 4
            first_matrix_element_index = stored_segment_end_offsets(column_index - 1) + 1
            diagonal_element_index = stored_segment_end_offsets(column_index)
            do matrix_element_index = first_matrix_element_index, diagonal_element_index
                row_index = matrix_partner_indices(matrix_element_index)
                if (first_column_debug_unit == -1) then
                    open(newunit=first_column_debug_unit, file="fort.96", status="unknown", form="formatted")
                end if
                write (first_column_debug_unit, *) "H(", row_index, column_index, ")= ", hamiltonian_values(matrix_element_index)
            end do
            ! ibeg = idiag + 1
        end do

        first_matrix_element_index = stored_segment_end_offsets(1) + 1
        do column_index = 2, block_csf_count, 4
            first_matrix_element_index = stored_segment_end_offsets(column_index - 1) + 1
            diagonal_element_index = stored_segment_end_offsets(column_index)
            do matrix_element_index = first_matrix_element_index, diagonal_element_index
                row_index = matrix_partner_indices(matrix_element_index)
                if (second_column_debug_unit == -1) then
                    open(newunit=second_column_debug_unit, file="fort.97", status="unknown", form="formatted")
                end if
                write (second_column_debug_unit, *) "H(", row_index, column_index, ")= ", hamiltonian_values(matrix_element_index)
            end do
            ! ibeg = idiag + 1
        end do

        first_matrix_element_index = stored_segment_end_offsets(2) + 1
        do column_index = 3, block_csf_count, 4
            first_matrix_element_index = stored_segment_end_offsets(column_index - 1) + 1
            diagonal_element_index = stored_segment_end_offsets(column_index)
            do matrix_element_index = first_matrix_element_index, diagonal_element_index
                row_index = matrix_partner_indices(matrix_element_index)
                if (third_column_debug_unit == -1) then
                    open(newunit=third_column_debug_unit, file="fort.98", status="unknown", form="formatted")
                end if
                write (third_column_debug_unit, *) "H(", row_index, column_index, ")= ", hamiltonian_values(matrix_element_index)
            end do
            ! ibeg = idiag + 1
        end do

        first_matrix_element_index = stored_segment_end_offsets(3) + 1
        do column_index = 4, block_csf_count, 4
            first_matrix_element_index = stored_segment_end_offsets(column_index - 1) + 1
            diagonal_element_index = stored_segment_end_offsets(column_index)
            do matrix_element_index = first_matrix_element_index, diagonal_element_index
                row_index = matrix_partner_indices(matrix_element_index)
                call ensure_debug_output()
                write (angular_debug_unit, *) "H(", row_index, column_index, ")= ", hamiltonian_values(matrix_element_index)
            end do
            ! ibeg = idiag + 1
        end do
        return

    end subroutine write_hamiltonian_records
end module write_hamiltonian_records_interface
