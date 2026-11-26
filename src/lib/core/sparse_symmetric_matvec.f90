! Multiply a symmetric Hamiltonian stored by sparse columns by a group of vectors.
module sparse_symmetric_matvec
    implicit none (type, external)

    ! Module interface
    private
    public :: multiply_sparse_symmetric_matrix
contains
    subroutine multiply_sparse_symmetric_matrix(matrix_order, vector_count, input_vectors, product_vectors)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use sparse_hamiltonian_storage, only: hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices
        use sparse_dot_and_scatter, only: sparse_dot_scatter
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: matrix_order
        integer, intent(in) :: vector_count
        real(real64), intent(in) :: input_vectors(matrix_order, vector_count)
        real(real64), intent(inout) :: product_vectors(matrix_order, vector_count)

        ! Local declarations
        integer :: column_index, vector_index
        integer(int64) :: column_start, column_end, column_entry_count
        real(real64) :: diagonal_plus_previous_terms, off_diagonal_dot_product
        product_vectors = 0.0_real64

        column_start = 1
        do column_index = 1, matrix_order
            ! The column's last stored entry is the diagonal. scatter the remaining entries into the other
            ! triangle.
            column_end = stored_segment_end_offsets(column_index)
            column_entry_count = column_end - column_start + 1
            do vector_index = 1, vector_count
                diagonal_plus_previous_terms = product_vectors(column_index, &
                    & vector_index) + hamiltonian_values(column_end)*input_vectors(column_index,vector_index)
                call sparse_dot_scatter (column_entry_count-1, input_vectors(:matrix_order,vector_index), &
                    & product_vectors(:matrix_order, &
                    & vector_index), matrix_partner_indices(column_start:column_end-1), &
                    & hamiltonian_values(column_start:column_end-1), input_vectors(&
                    column_index,vector_index), off_diagonal_dot_product)
                product_vectors(column_index,vector_index) = diagonal_plus_previous_terms + off_diagonal_dot_product
            end do
            column_start = column_end + 1
        end do

        return
    end subroutine multiply_sparse_symmetric_matrix
end module sparse_symmetric_matvec
