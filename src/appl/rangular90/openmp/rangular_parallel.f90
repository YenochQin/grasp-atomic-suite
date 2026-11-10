! Bounded batches of independent CSF pairs; only the initial thread commits
! results to the original MCP writer. Worker module state is TLS, never copied.
module rangular_parallel
   use iso_fortran_env, only: int64, real64, error_unit
   use omp_lib, only: omp_in_parallel, omp_get_max_threads, omp_get_thread_limit, omp_set_dynamic, &
                     omp_get_wtime, omp_get_num_threads, omp_get_thread_num
   use vast_kind_param, only: double
   use parameter_def, only: orbital_capacity => NNNW
   use buffer_C, only: raw_labels => LABEL, raw_coefficients => COEFF, &
                       raw_coefficient_count => NVCOEF, pair_buffer_capacity => NBDIM
   use orb_C, only: current_csf_count => NCF
   use mcp_C, only: diagonal_only => DIAG, first_order_only => LFORDR
   use iccu_C, only: reference_limits => ICCUT
   use debug_C, only: angular_debug_flags => LDBPA, coefficient_debug => IBUG1
   use onescalar_I, only: onescalar
   use rkco_gg_I, only: rkco_gg
   use cord_I, only: cord
   use alcbuf_I, only: alcbuf
   use rangular_profile, only: record_seconds_stage
   implicit none
   private
   public :: parallel_pairs_available, configure_pair_engine, begin_pair_block, end_pair_block, load_parallel_pair
   logical, parameter :: parallel_pairs_available = .true.
   integer, parameter :: bytes_per_coefficient = 5*storage_size(0)/8 + storage_size(0.0_double)/8

   type :: pair_tile
      integer :: first_row = 0, last_row = -1, stored_coefficients = 0, scratch_unit = 0
      logical :: spilled = .false.
      character(len=256) :: error_message = ''
      integer, allocatable :: first_orbitals(:), second_orbitals(:), coefficient_counts(:), coefficient_offsets(:)
      integer(int64), allocatable :: stream_positions(:)
      real(double), allocatable :: one_electron_coefficients(:), coefficients(:)
      integer, allocatable :: labels(:,:)
   end type pair_tile

   type(pair_tile), allocatable, save :: result_tiles(:)
   integer(int64), allocatable, save :: worker_chunk_counts(:)
   integer, save :: worker_count = 1, requested_chunk_size = 256, result_memory_mib = 64
   integer, save :: reference_limit = 0, cached_column = 0, first_cached_row = 0, last_cached_row = -1
   integer, save :: rows_per_tile = 1, used_tiles = 0
   integer(int64), save :: block_pair_count = 0, block_spill_count = 0
   real(real64), save :: batch_compute_seconds = 0.0_real64
contains
   subroutine configure_pair_engine(requested_threads, pair_chunk_size, memory_limit_mib)
      integer, intent(in) :: requested_threads, pair_chunk_size, memory_limit_mib
      integer :: maximum_worker_count
      if (omp_in_parallel()) error stop 'rangular pair engine cannot be nested in another parallel region'
      maximum_worker_count = omp_get_thread_limit()
      worker_count = min(omp_get_max_threads(), maximum_worker_count)
      if (requested_threads > 0) worker_count = requested_threads
      if (worker_count < 1 .or. worker_count > maximum_worker_count) error stop 'Invalid OpenMP thread count/limit'
      if (int(worker_count,int64)*4_int64 > int(huge(worker_count),int64)) error stop 'Thread count is too large'
      if (pair_chunk_size < 1 .or. memory_limit_mib < 0) error stop 'Invalid parallel buffer/chunk option'
      requested_chunk_size = pair_chunk_size
      result_memory_mib = memory_limit_mib
      ! Stable teams/affinity, no nested regions/tasks/teams or resource pauses:
      ! OpenMP 5.2 section 5.2 guarantees TLS/cache/pointer persistence between
      ! these regions. Buffers are initialized and released once per block.
      call omp_set_dynamic(.false.)
      write(*,'(A,I0,A,I0,A,I0)') 'OpenMP pair engine: threads=', worker_count, &
         ', maximum pairs/chunk=', requested_chunk_size, ', result buffer MiB=', result_memory_mib
   end subroutine configure_pair_engine

   subroutine begin_pair_block(block_index)
      integer, intent(in) :: block_index
      integer :: tile_index, tile_count, coefficient_capacity, allocation_status
      integer(int64) :: capacity_wide
      if (allocated(result_tiles)) error stop 'Previous OpenMP block was not released'
      if (any(angular_debug_flags) .or. coefficient_debug /= 0) error stop 'Use serial rangular for angular debug output'
      reference_limit = reference_limits(block_index)
      tile_count = min(4*worker_count, current_csf_count)
      if (tile_count < 1) error stop 'Empty OpenMP block'
      capacity_wide = int(result_memory_mib,int64)*1048576_int64/(int(tile_count,int64)*bytes_per_coefficient)
      if (capacity_wide > int(huge(coefficient_capacity),int64)) error stop 'Result buffer capacity exceeds integer range'
      coefficient_capacity = int(capacity_wide)
      allocate(result_tiles(tile_count), worker_chunk_counts(worker_count), stat=allocation_status)
      if (allocation_status /= 0) error stop 'Cannot allocate OpenMP result controller'
      worker_chunk_counts = 0
      do tile_index = 1, tile_count
         associate(tile => result_tiles(tile_index))
            allocate(tile%first_orbitals(min(requested_chunk_size,current_csf_count)), &
                     tile%second_orbitals(min(requested_chunk_size,current_csf_count)), &
                     tile%one_electron_coefficients(min(requested_chunk_size,current_csf_count)), &
                     tile%coefficient_counts(min(requested_chunk_size,current_csf_count)), &
                     tile%coefficient_offsets(min(requested_chunk_size,current_csf_count)), &
                     tile%stream_positions(min(requested_chunk_size,current_csf_count)), &
                     tile%labels(5,coefficient_capacity), tile%coefficients(coefficient_capacity), stat=allocation_status)
            if (allocation_status /= 0) error stop 'Cannot allocate bounded OpenMP result buffers'
         end associate
      end do
      cached_column = 0
      first_cached_row = 0
      last_cached_row = -1
      used_tiles = 0
      block_pair_count = 0_int64
      block_spill_count = 0_int64
      batch_compute_seconds = 0.0_real64
      !$omp parallel num_threads(worker_count) default(none) shared(worker_count)
      if (omp_get_num_threads() /= worker_count) error stop 'OpenMP runtime did not provide the requested fixed team'
      ! MCP already allocated the initial thread's buffer.
      if (omp_get_thread_num() /= 0) call alcbuf(1)
      raw_coefficient_count = 0
      !$omp end parallel
   end subroutine begin_pair_block

   subroutine clear_tile_files()
      integer :: tile_index, close_status
      if (.not. allocated(result_tiles)) return
      do tile_index = 1, size(result_tiles)
         if (result_tiles(tile_index)%spilled) then
            close(result_tiles(tile_index)%scratch_unit, status='delete', iostat=close_status)
            if (close_status /= 0) error stop 'Cannot close OpenMP scratch result stream'
         end if
         result_tiles(tile_index)%spilled = .false.
         result_tiles(tile_index)%error_message = ''
         result_tiles(tile_index)%stored_coefficients = 0
      end do
   end subroutine clear_tile_files

   subroutine spill_previous_pairs(tile, previous_pair_count)
      type(pair_tile), intent(inout) :: tile
      integer, intent(in) :: previous_pair_count
      integer :: pair_index, first_coefficient, last_coefficient, io_status
      open(newunit=tile%scratch_unit, status='scratch', access='stream', form='unformatted', &
           action='readwrite', iostat=io_status)
      if (io_status /= 0) then
         tile%error_message = 'Cannot open scratch stream for bounded pair results'
         return
      end if
      tile%spilled = .true.
      do pair_index = 1, previous_pair_count
         if (tile%coefficient_counts(pair_index) == 0) cycle
         first_coefficient = tile%coefficient_offsets(pair_index)
         last_coefficient = first_coefficient + tile%coefficient_counts(pair_index) - 1
         write(tile%scratch_unit, iostat=io_status) tile%labels(:,first_coefficient:last_coefficient), &
                                                  tile%coefficients(first_coefficient:last_coefficient)
         if (io_status /= 0) then
            tile%error_message = 'Cannot spill earlier pair results'
            return
         end if
      end do
   end subroutine spill_previous_pairs

   subroutine compute_tile(column_index, tile)
      integer, intent(in) :: column_index
      type(pair_tile), intent(inout) :: tile
      integer :: row_index, pair_index, first_coefficient, last_coefficient, io_status
      integer(int64) :: stream_position
      real(double) :: shell_coefficients(orbital_capacity)
      stream_position = 1_int64
      do row_index = tile%first_row, tile%last_row
         pair_index = row_index - tile%first_row + 1
         tile%first_orbitals(pair_index) = 0
         tile%second_orbitals(pair_index) = 0
         tile%one_electron_coefficients(pair_index) = 0.0_double
         tile%coefficient_counts(pair_index) = 0
         tile%coefficient_offsets(pair_index) = tile%stored_coefficients + 1
         tile%stream_positions(pair_index) = stream_position
         if (diagonal_only .or. (first_order_only .and. row_index > reference_limit)) then
            if (row_index /= column_index) cycle
         end if
         if (row_index /= column_index) then
            call onescalar(column_index, row_index, tile%first_orbitals(pair_index), &
                           tile%second_orbitals(pair_index), shell_coefficients)
            ! ONESCALAR's early rejection defines only its first orbital tag.
            if (tile%first_orbitals(pair_index) /= 0) then
               tile%one_electron_coefficients(pair_index) = shell_coefficients(1)
            else
               tile%second_orbitals(pair_index) = 0
            end if
         end if
         raw_coefficient_count = 0
         call rkco_gg(column_index, row_index, cord, 0, 1)
         tile%coefficient_counts(pair_index) = raw_coefficient_count
         if (raw_coefficient_count == 0) cycle
         if (.not. tile%spilled) then
            if (raw_coefficient_count > size(tile%coefficients)-tile%stored_coefficients) &
               call spill_previous_pairs(tile, pair_index-1)
            if (len_trim(tile%error_message) /= 0) return
         end if
         if (tile%spilled) then
            write(tile%scratch_unit, iostat=io_status) raw_labels(1:5,:raw_coefficient_count), &
                                                     raw_coefficients(:raw_coefficient_count)
            if (io_status /= 0) then
               tile%error_message = 'Cannot write current pair results'
               return
            end if
         else
            first_coefficient = tile%stored_coefficients + 1
            last_coefficient = tile%stored_coefficients + raw_coefficient_count
            tile%labels(:,first_coefficient:last_coefficient) = raw_labels(1:5,:raw_coefficient_count)
            tile%coefficients(first_coefficient:last_coefficient) = raw_coefficients(:raw_coefficient_count)
            tile%stored_coefficients = last_coefficient
         end if
         stream_position = stream_position + int(raw_coefficient_count,int64)*bytes_per_coefficient
      end do
   end subroutine compute_tile

   subroutine prepare_batch(column_index, row_index)
      integer, intent(in) :: column_index, row_index
      integer :: tile_index, batch_rows, permitted_last_row, thread_index
      real(real64) :: start_seconds
      call clear_tile_files()
      cached_column = column_index
      first_cached_row = row_index
      permitted_last_row = column_index
      if (first_order_only .and. row_index <= reference_limit) permitted_last_row = min(column_index,reference_limit)
      rows_per_tile = min(requested_chunk_size, 1+(permitted_last_row-row_index)/size(result_tiles))
      batch_rows = int(min(int(permitted_last_row-row_index+1,int64), &
                          int(size(result_tiles),int64)*rows_per_tile))
      last_cached_row = row_index + batch_rows - 1
      used_tiles = 1 + (batch_rows-1)/rows_per_tile
      do tile_index = 1, used_tiles
         result_tiles(tile_index)%first_row = row_index + (tile_index-1)*rows_per_tile
         result_tiles(tile_index)%last_row = min(last_cached_row, result_tiles(tile_index)%first_row+rows_per_tile-1)
      end do
      start_seconds = omp_get_wtime()
      !$omp parallel num_threads(worker_count) default(none) &
      !$omp shared(result_tiles,used_tiles,column_index,worker_count,worker_chunk_counts) private(tile_index,thread_index)
      if (omp_get_num_threads() /= worker_count) error stop 'OpenMP runtime did not provide the requested fixed team'
      thread_index = omp_get_thread_num() + 1
      !$omp do schedule(dynamic,1)
      do tile_index = 1, used_tiles
         call compute_tile(column_index,result_tiles(tile_index))
         worker_chunk_counts(thread_index) = worker_chunk_counts(thread_index) + 1_int64
      end do
      !$omp end do
      !$omp end parallel
      batch_compute_seconds = batch_compute_seconds + omp_get_wtime() - start_seconds
      block_pair_count = block_pair_count + int(batch_rows,int64)
      do tile_index = 1, used_tiles
         if (len_trim(result_tiles(tile_index)%error_message) /= 0) then
            write(error_unit,'(A,I0,A,I0,2A)') 'OpenMP pair result failure at column ',column_index, &
               ', first row ',result_tiles(tile_index)%first_row, ': ',trim(result_tiles(tile_index)%error_message)
            error stop 1
         end if
         if (result_tiles(tile_index)%spilled) block_spill_count = block_spill_count + 1_int64
      end do
   end subroutine prepare_batch

   subroutine load_parallel_pair(column_index, row_index, first_orbital, second_orbital, one_electron_coefficient)
      integer, intent(in) :: column_index, row_index
      integer, intent(out) :: first_orbital, second_orbital
      real(double), intent(out) :: one_electron_coefficient
      integer :: tile_index, pair_index, first_coefficient, last_coefficient, io_status
      if (.not. allocated(result_tiles)) error stop 'OpenMP pair block has not been initialized'
      if (cached_column /= column_index .or. row_index < first_cached_row .or. row_index > last_cached_row) &
         call prepare_batch(column_index,row_index)
      tile_index = 1 + (row_index-first_cached_row)/rows_per_tile
      associate(tile => result_tiles(tile_index))
         pair_index = row_index-tile%first_row+1
         first_orbital = tile%first_orbitals(pair_index)
         second_orbital = tile%second_orbitals(pair_index)
         one_electron_coefficient = tile%one_electron_coefficients(pair_index)
         raw_coefficient_count = tile%coefficient_counts(pair_index)
         do while (pair_buffer_capacity < raw_coefficient_count)
            call alcbuf(2)
         end do
         if (raw_coefficient_count == 0) return
         if (tile%spilled) then
            read(tile%scratch_unit,pos=tile%stream_positions(pair_index),iostat=io_status) &
               raw_labels(1:5,:raw_coefficient_count),raw_coefficients(:raw_coefficient_count)
            if (io_status /= 0) error stop 'Cannot read back OpenMP pair result stream'
         else
            first_coefficient = tile%coefficient_offsets(pair_index)
            last_coefficient = first_coefficient+raw_coefficient_count-1
            raw_labels(1:5,:raw_coefficient_count) = tile%labels(:,first_coefficient:last_coefficient)
            raw_coefficients(:raw_coefficient_count) = tile%coefficients(first_coefficient:last_coefficient)
         end if
      end associate
   end subroutine load_parallel_pair

   subroutine end_pair_block(block_index)
      integer, intent(in) :: block_index
      character(len=24) :: result_storage
      result_storage = 'omp-memory'
      if (block_spill_count > 0) result_storage = 'omp-spill'
      call record_seconds_stage('parallel_batch_compute',block_index,batch_compute_seconds,block_pair_count,result_storage)
      write(*,'(A,I0,A,I0,A,I0)') 'OpenMP block ',block_index,': active workers=', &
         count(worker_chunk_counts > 0),', spilled chunks=',block_spill_count
      write(*,'(A,*(I0,1X))') 'OpenMP chunks per worker: ',worker_chunk_counts
      call clear_tile_files()
      deallocate(result_tiles,worker_chunk_counts)
      !$omp parallel num_threads(worker_count) default(none) shared(worker_count)
      if (omp_get_num_threads() /= worker_count) error stop 'OpenMP team changed during block cleanup'
      ! MCP releases the initial thread's buffer after this routine returns.
      if (omp_get_thread_num() /= 0) call alcbuf(3)
      !$omp end parallel
   end subroutine end_pair_block
end module rangular_parallel
