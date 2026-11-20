module mpi_sparse_matrix_product
    implicit none (type, external)
    private
    public :: spicmvmpi
contains
      ! Apply the symmetric sparse Hamiltonian to a batch of Davidson vectors.
      subroutine spicmvmpi (matrix_order,vector_count,input_vectors,product_vectors)

      use, intrinsic :: iso_fortran_env, only: real64
      use mpi_state, only: mpi_process_count, mpi_rank
      use sparse_hamiltonian_storage, only: hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices
      use real_vector_initialization, only: dinit
      use sparse_dot_and_scatter, only: dmerge
      use mpi_utilities, only: gdsummpi
      implicit none (type, external)
      integer, intent(in) :: matrix_order
      integer, intent(in) :: vector_count
      real(real64), intent(in) :: input_vectors(matrix_order,vector_count)
      real(real64), intent(inout) :: product_vectors(matrix_order,vector_count)
      integer :: column_index, segment_start, segment_end, segment_element_count, vector_index
      real(real64) :: partial_product, off_diagonal_product
      call dinit (matrix_order*vector_count, 0.0_real64, product_vectors, 1)

      ! Ranks own columns round-robin. Each stored segment ends with its diagonal entry.
      segment_start = 1
      do column_index = mpi_rank + 1, matrix_order, mpi_process_count
         segment_end = stored_segment_end_offsets(column_index)
         segment_element_count = segment_end - segment_start + 1
         do vector_index = 1, vector_count
            partial_product = product_vectors(column_index,vector_index) + &
                              hamiltonian_values(segment_end)*input_vectors(column_index,vector_index)
            ! Accumulate both symmetric contributions from each stored off-diagonal entry.
            call dmerge (segment_element_count-1, input_vectors(:matrix_order,vector_index), &
                         product_vectors(:matrix_order,vector_index), matrix_partner_indices(segment_start:segment_end-1), &
                         hamiltonian_values(segment_start:segment_end-1), input_vectors(column_index,vector_index), &
                         off_diagonal_product)
            product_vectors(column_index,vector_index) = partial_product + off_diagonal_product
         end do
         segment_start = segment_end + 1
      end do

      ! Sum the rank-local products and return the full result on every rank.
      call gdsummpi (product_vectors, matrix_order*vector_count)

      return
      end subroutine spicmvmpi
end module mpi_sparse_matrix_product
