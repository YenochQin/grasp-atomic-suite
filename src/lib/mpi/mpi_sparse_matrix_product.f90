module mpi_sparse_matrix_product
    implicit none (type, external)

    ! Module interface
    private
    public :: multiply_mpi_sparse_matrix
contains
    ! Apply the symmetric sparse Hamiltonian to a batch of Davidson vectors.
    subroutine multiply_mpi_sparse_matrix (matrix_order,vector_count,input_vectors,product_vectors)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use mpi_state, only: mpi_process_count, mpi_rank
        use sparse_hamiltonian_storage, only: hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices
        use real_vector_initialization, only: initialize_real_vector
        use sparse_dot_and_scatter, only: sparse_dot_scatter
        use mpi_f08, only: mpi_count_kind
        use mpi_element_counts, only: checked_mpi_count
        use mpi_utilities, only: allreduce_mpi_real_sum
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: matrix_order
        integer, intent(in) :: vector_count
        real(real64), intent(in) :: input_vectors(matrix_order,vector_count)
        real(real64), intent(inout) :: product_vectors(matrix_order,vector_count)

        ! Local declarations
        integer :: value_count
        integer :: column_index, segment_start, segment_end, segment_element_count, vector_index
        real(real64) :: partial_product, off_diagonal_product
        value_count = checked_mpi_count(size(product_vectors,kind=mpi_count_kind), "multiply_mpi_sparse_matrix")
        call initialize_real_vector (value_count, 0.0_real64, product_vectors, 1)

        ! Ranks own columns round-robin. Each stored segment ends with its diagonal entry.
        ! The offset array is indexed by global column. Consecutive local segments are adjacent
        ! even when their global column indices differ by more than one.
        segment_start = 1
        do column_index = mpi_rank + 1, matrix_order, mpi_process_count
            segment_end = stored_segment_end_offsets(column_index)
            segment_element_count = segment_end - segment_start + 1
            do vector_index = 1, vector_count
                partial_product = product_vectors(column_index,vector_index) + &
                    hamiltonian_values(segment_end)*input_vectors(column_index,vector_index)
                ! Accumulate both symmetric contributions from each stored off-diagonal entry.
                call sparse_dot_scatter (segment_element_count-1, input_vectors(:matrix_order,vector_index), &
                    product_vectors(:matrix_order,vector_index), matrix_partner_indices(segment_start:segment_end-1), &
                    hamiltonian_values(segment_start:segment_end-1), input_vectors(column_index,vector_index), &
                    off_diagonal_product)
                product_vectors(column_index,vector_index) = partial_product + off_diagonal_product
            end do
            segment_start = segment_end + 1
        end do

        ! Sum the rank-local products and return the full result on every rank.
        call allreduce_mpi_real_sum (product_vectors, value_count)

        return
    end subroutine multiply_mpi_sparse_matrix
end module mpi_sparse_matrix_product
