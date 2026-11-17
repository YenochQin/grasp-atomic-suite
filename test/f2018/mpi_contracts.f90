program mpi_contracts
    use mpi_c, only: host => processor_name, ierr => mpi_error_code, istat => mpi_status_buffer, &
       & lenhost => processor_name_length, myid => mpi_rank, nprocs => mpi_process_count, mpi_abort, mpi_allreduce, &
       & mpi_bcast, mpi_character, mpi_comm_world, mpi_double_precision, mpi_finalize, mpi_logical, mpi_success, &
       & mpi_sum, mpi_type_free, mpi_type_match_size, mpi_type_size, mpi_typeclass_integer, mpix_bytes, startmpi
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use, intrinsic :: iso_fortran_env, only: int8
    implicit none (type, external)
    integer :: byte_type, wide_type, width
    integer(int8) :: packed(2,3)
    integer(int64) :: wide(2)
    real(real64) :: values(2), local, total
    logical :: enabled
    character(12) :: message

    call startmpi(myid,nprocs,host,lenhost)
    if (lenhost <= 0 .or. lenhost > len(host)) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    call mpix_bytes(1,byte_type,ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,2,ierr)
    call MPI_Type_size(byte_type,width,ierr)
    if (width /= 1) call MPI_Abort(MPI_COMM_WORLD,3,ierr)
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(wide(1))/8,wide_type,ierr)
    packed = 0
    wide = 0
    values = 0
    enabled = .false.
    message = ""
    if (myid == 0) then
        packed = reshape([1_int8,2_int8,3_int8,4_int8,5_int8,6_int8],[2,3])
        wide = [2_int64**40,2_int64**40+7_int64]
        values = [0.25_real64,-1.5_real64]
        enabled = .true.
        message = "rank-zero"
    end if
    call MPI_Bcast(packed,size(packed),byte_type,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(wide,size(wide),wide_type,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(values,size(values),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(enabled,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(message,len(message),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    if (any(packed /= reshape([1_int8,2_int8,3_int8,4_int8,5_int8,6_int8],[2,3]))) &
        call MPI_Abort(MPI_COMM_WORLD,4,ierr)
    if (any(wide /= [2_int64**40,2_int64**40+7_int64])) call MPI_Abort(MPI_COMM_WORLD,5,ierr)
    if (any(values /= [0.25_real64,-1.5_real64])) call MPI_Abort(MPI_COMM_WORLD,6,ierr)
    if (.not.enabled .or. message /= "rank-zero") call MPI_Abort(MPI_COMM_WORLD,7,ierr)
    local = real(myid+1,real64)
    call MPI_Allreduce(local,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (total /= real(nprocs*(nprocs+1)/2,real64)) call MPI_Abort(MPI_COMM_WORLD,8,ierr)
    call MPI_Type_free(byte_type,ierr)
    call MPI_Finalize(ierr)
end program mpi_contracts
