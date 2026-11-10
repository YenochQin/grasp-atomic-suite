! Same controller interface as rangular_omp; the original serial kernel remains the baseline.
module rangular_parallel
   use vast_kind_param, only: double
   implicit none
   private
   public :: parallel_pairs_available, configure_pair_engine, begin_pair_block, end_pair_block, load_parallel_pair
   logical, parameter :: parallel_pairs_available = .false.
contains
   subroutine configure_pair_engine(requested_threads, pair_chunk_size, result_memory_mib)
      integer, intent(in) :: requested_threads, pair_chunk_size, result_memory_mib
      if (requested_threads > 1 .or. pair_chunk_size /= 256 .or. result_memory_mib /= 64) &
         error stop 'Parallel options require the rangular_omp executable'
   end subroutine configure_pair_engine

   subroutine begin_pair_block(block_index)
      integer, intent(in) :: block_index
      if (block_index < 1) error stop 'Invalid block index'
   end subroutine begin_pair_block

   subroutine end_pair_block(block_index)
      integer, intent(in) :: block_index
      if (block_index < 1) error stop 'Invalid block index'
   end subroutine end_pair_block

   subroutine load_parallel_pair(column_index, row_index, first_orbital, second_orbital, one_electron_coefficient)
      integer, intent(in) :: column_index, row_index
      integer, intent(out) :: first_orbital, second_orbital
      real(double), intent(out) :: one_electron_coefficient
      first_orbital = 0
      second_orbital = 0
      one_electron_coefficient = 0.0_double
      if (column_index >= row_index) error stop 'Serial target must use the original pair kernel'
      error stop 'Invalid CSF pair'
   end subroutine load_parallel_pair
end module rangular_parallel
