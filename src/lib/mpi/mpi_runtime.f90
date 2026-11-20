module mpi_runtime
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_mpi_runtime
    public :: mpix_bytes
    public :: finalize_mpi_runtime
contains
    ! Initialize MPI_COMM_WORLD and return the rank, process count and processor name.
    subroutine initialize_mpi_runtime (rank, process_count, processor_name, processor_name_length)

        ! Module dependencies
        use mpi, only: mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_get_processor_name, mpi_init, mpi_success
        implicit none (type, external)

        ! Arguments
        integer, intent(out)  :: rank,  process_count, processor_name_length
        character(len=*), intent(out) :: processor_name

        ! Local declarations
        integer :: mpi_return_code

        call mpi_init (mpi_return_code)
        call mpi_comm_rank (mpi_comm_world, rank, mpi_return_code)
        call mpi_comm_size (mpi_comm_world, process_count, mpi_return_code)
        call mpi_get_processor_name (processor_name, processor_name_length, mpi_return_code)
        if (mpi_return_code /= mpi_success) print *, " ierr=", mpi_return_code, " id=", rank

        return
    end subroutine initialize_mpi_runtime
    ! Create and commit a datatype consisting of byte_count contiguous MPI_BYTE entries.
    subroutine mpix_bytes (byte_count, byte_datatype, mpi_return_code)

        ! Module dependencies
        use mpi, only: mpi_byte, mpi_type_commit, mpi_type_contiguous
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: byte_count
        integer, intent(out) :: byte_datatype, mpi_return_code

        ! Local declarations
        integer :: create_status

        call mpi_type_contiguous (byte_count, mpi_byte, byte_datatype, create_status)
        call mpi_type_commit (byte_datatype, mpi_return_code)

        mpi_return_code = mpi_return_code + create_status

        return
    end subroutine mpix_bytes
    ! Print the caller's message, synchronize MPI_COMM_WORLD and finalize MPI.
    ! All ranks must reach this routine to complete the barrier before finalizing.
    subroutine finalize_mpi_runtime (message, rank)

        ! Module dependencies
        use mpi, only: mpi_barrier, mpi_comm_world, mpi_finalize
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  ::  rank
        character(len=*), intent(in) :: message

        ! Local declarations
        integer :: mpi_return_code

        print *, "mpi stopped by node-", rank, " from ", message
        call mpi_barrier (mpi_comm_world,mpi_return_code)
        call mpi_finalize (mpi_return_code)

        stop
    end subroutine finalize_mpi_runtime
end module mpi_runtime
