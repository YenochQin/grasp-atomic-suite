program mpi_contracts
    use mpi_state, only: processor_name, mpi_error_code, processor_name_length, mpi_rank, &
        mpi_process_count, mpi_abort, mpi_allreduce, mpi_bcast, mpi_character, mpi_comm_world, mpi_double_precision, &
        mpi_finalize, mpi_logical, mpi_success, mpi_sum, mpi_type_free, mpi_type_match_size, mpi_type_size, &
        mpi_count_kind, mpi_datatype, mpi_get_count, mpi_integer, mpi_sendrecv, mpi_status, mpi_typeclass_integer, &
        create_contiguous_mpi_byte_datatype, initialize_mpi_runtime
    use mpi_f08, only: mpi_allgather, mpi_allgatherv
    use mpi_element_counts, only: checked_mpi_count, max_mpi_count, mpi_gather_displacements
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use, intrinsic :: iso_fortran_env, only: int8
    implicit none (type, external)
    type(mpi_datatype) :: byte_type, wide_type
    type(mpi_status) :: receive_status
    integer :: width, received_rank, previous_rank, next_rank, local_count, total_count, rank_index
    integer, allocatable :: counts(:), offsets(:), gathered(:)
    integer(int8) :: packed(2,3)
    integer(int64) :: wide(2)
    real(real64) :: values(2), local, total
    logical :: enabled
    character(12) :: message, scenario

    call initialize_mpi_runtime(mpi_rank,mpi_process_count,processor_name,processor_name_length)
    call get_command_argument(1, scenario)
    if (scenario == "overflow" .and. mpi_rank == 1) then
        local_count = checked_mpi_count(max_mpi_count + 1_mpi_count_kind, "MPI count probe on rank 1")
        print *, "Unexpectedly accepted count ", local_count
    end if
    if (checked_mpi_count(max_mpi_count, "MPI boundary probe") /= max_mpi_count) &
        call MPI_Abort(MPI_COMM_WORLD,9,mpi_error_code)
    previous_rank = modulo(mpi_rank-1, mpi_process_count)
    next_rank = modulo(mpi_rank+1, mpi_process_count)
    call mpi_sendrecv(mpi_rank, 1, mpi_integer, next_rank, 0, received_rank, 1, mpi_integer, previous_rank, 0, &
        mpi_comm_world, receive_status, mpi_error_code)
    call mpi_get_count(receive_status, mpi_integer, width, mpi_error_code)
    if (width /= 1 .or. received_rank /= previous_rank .or. receive_status%mpi_source /= previous_rank) &
        call MPI_Abort(MPI_COMM_WORLD,10,mpi_error_code)

    ! Match the RMCDHF variable-count exchange, including zero-length contributions.
    local_count = merge(1, 0, modulo(mpi_rank,2) == 0)
    allocate(counts(mpi_process_count), offsets(mpi_process_count))
    call mpi_allgather(local_count, 1, mpi_integer, counts, 1, mpi_integer, mpi_comm_world, mpi_error_code)
    call mpi_gather_displacements(counts, offsets, total_count)
    allocate(gathered(total_count))
    call mpi_allgatherv([mpi_rank], local_count, mpi_integer, gathered, counts, offsets, mpi_integer, &
        mpi_comm_world, mpi_error_code)
    if (total_count /= (mpi_process_count+1)/2 .or. &
        any(gathered /= [(rank_index, rank_index=0,mpi_process_count-1,2)])) call MPI_Abort(MPI_COMM_WORLD,11,mpi_error_code)
    if (processor_name_length <= 0 .or. processor_name_length > len(processor_name)) call MPI_Abort(MPI_COMM_WORLD,1,mpi_error_code)
    call create_contiguous_mpi_byte_datatype(1,byte_type,mpi_error_code)
    if (mpi_error_code /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,2,mpi_error_code)
    call MPI_Type_size(byte_type,width,mpi_error_code)
    if (width /= 1) call MPI_Abort(MPI_COMM_WORLD,3,mpi_error_code)
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(wide(1))/8,wide_type,mpi_error_code)
    packed = 0
    wide = 0
    values = 0
    enabled = .false.
    message = ""
    if (mpi_rank == 0) then
        packed = reshape([1_int8,2_int8,3_int8,4_int8,5_int8,6_int8],[2,3])
        wide = [2_int64**40,2_int64**40+7_int64]
        values = [0.25_real64,-1.5_real64]
        enabled = .true.
        message = "rank-zero"
    end if
    call MPI_Bcast(packed,size(packed),byte_type,0,MPI_COMM_WORLD,mpi_error_code)
    call MPI_Bcast(wide,size(wide),wide_type,0,MPI_COMM_WORLD,mpi_error_code)
    call MPI_Bcast(values,size(values),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_error_code)
    call MPI_Bcast(enabled,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_error_code)
    call MPI_Bcast(message,len(message),MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_error_code)
    if (any(packed /= reshape([1_int8,2_int8,3_int8,4_int8,5_int8,6_int8],[2,3]))) &
        call MPI_Abort(MPI_COMM_WORLD,4,mpi_error_code)
    if (any(wide /= [2_int64**40,2_int64**40+7_int64])) call MPI_Abort(MPI_COMM_WORLD,5,mpi_error_code)
    if (any(values /= [0.25_real64,-1.5_real64])) call MPI_Abort(MPI_COMM_WORLD,6,mpi_error_code)
    if (.not.enabled .or. message /= "rank-zero") call MPI_Abort(MPI_COMM_WORLD,7,mpi_error_code)
    local = real(mpi_rank+1,real64)
    call MPI_Allreduce(local,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_error_code)
    if (total /= real(mpi_process_count*(mpi_process_count+1)/2,real64)) call MPI_Abort(MPI_COMM_WORLD,8,mpi_error_code)
    call MPI_Type_free(byte_type,mpi_error_code)
    call MPI_Finalize(mpi_error_code)
end program mpi_contracts
