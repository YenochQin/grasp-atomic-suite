module mpi_utilities
    ! Module dependencies
    use mpi_runtime, only: initialize_mpi_runtime, create_contiguous_mpi_byte_datatype, finalize_mpi_runtime
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_mpi_runtime
    public :: start_mpi_application
    public :: initialize_mpi_runtime_with_status
    public :: finish_mpi_application
    public :: finalize_mpi_runtime
    public :: print_mpi_message
    public :: check_mpi_checkpoint
    public :: create_contiguous_mpi_byte_datatype
    public :: allreduce_mpi_real_sum
    public :: reduce_mpi_real_sum_to_root
    public :: reduce_mpi_real_max_to_root
    public :: allreduce_mpi_integer_sum
contains

    ! Start MPI in the caller's working directory and optionally redirect worker output.
    subroutine start_mpi_application (rank, process_count, processor_name, processor_name_length, start_ticks, &
        program_name, description, input_file_description, &
        output_file_description, quiet_workers)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: error_unit, int64, output_unit
        use mpi_run_logging, only: report_mpi_rank_values
        use run_timer, only: begin_run_timer

        implicit none (type, external)

        ! Arguments
        integer, intent(out) :: rank, process_count, processor_name_length
        integer(int64), intent(out) :: start_ticks
        character(len=*), intent(out) :: processor_name
        character(len=*), intent(in) :: program_name, description
        character(len=*), intent(in) :: input_file_description, output_file_description
        logical, intent(in) :: quiet_workers

        ! Local declarations
        integer :: output_open_status
        character(len=256) :: output_message

        call initialize_mpi_runtime (rank, process_count, processor_name, processor_name_length)
        start_ticks = 0_int64
        if (rank == 0) call begin_run_timer (start_ticks, program_name)

        if (rank == 0) then
            write (output_unit, *) trim(description)
            write (output_unit, *) "Input files:  " // trim(input_file_description)
            write (output_unit, *) "Output files: " // trim(output_file_description)
        else if (quiet_workers) then
            ! Reconnect the standard stream used by WRITE(*) and PRINT.
            open (unit=output_unit, file="/dev/null", status="OLD", action="WRITE", &
                iostat=output_open_status, iomsg=output_message)
            ! A failed redirect is a warning. The worker continues to participate in collectives.
            if (output_open_status /= 0) then
                write(error_unit,"(a,i0,a,i0,2a)") "MPI worker output warning: rank=", &
                    rank, "; iostat=", output_open_status, "; ", trim(output_message)
            end if
        end if

        call report_mpi_rank_values (rank, process_count, "Participating nodes:", processor_name, as_hosts=.true.)
        return
    end subroutine start_mpi_application

    ! Start MPI, report the participating ranks and set up rank-local working directories.
    subroutine initialize_mpi_runtime_with_status (rank, process_count, processor_name, processor_name_length, start_ticks, &
        startup_directory, serial_io_directory, rank_scratch_root, program_name)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: int64, output_unit
        use mpi_run_logging, only: report_mpi_rank_values
        use run_timer, only: begin_run_timer
        use mpi_work_paths, only: prepare_mpi_work_directories

        implicit none (type, external)

        ! Arguments
        integer, intent(out)  :: rank,  process_count, processor_name_length
        integer(int64), intent(out) :: start_ticks
        character(len=*), intent(out) :: processor_name
        character(len=*), intent(in) :: program_name
        character(len=*), intent(out) :: startup_directory, serial_io_directory, rank_scratch_root

        call initialize_mpi_runtime (rank, process_count, processor_name, processor_name_length)
        start_ticks = 0_int64
        if (rank == 0) call begin_run_timer (start_ticks, program_name)
        call report_mpi_rank_values (rank, process_count, "Participating nodes:", processor_name, as_hosts=.true.)

        ! Keep the serial I/O directory for rank 0 and select a scratch root for each rank.
        call prepare_mpi_work_directories (startup_directory, serial_io_directory, rank_scratch_root)

        call report_mpi_rank_values (rank, process_count, "Start Dir:", startup_directory, as_hosts=.false.)
        if (rank == 0) then
            write (output_unit, "(2a)") "Serial I/O Dir (rank 0): ", trim(serial_io_directory)
        end if
        call report_mpi_rank_values (rank, process_count, "Work Dir (Parallel I/O):", rank_scratch_root, as_hosts=.false.)
        return
    end subroutine initialize_mpi_runtime_with_status

    ! All ranks finalize MPI before rank 0 reports the end of the run.
    subroutine finish_mpi_application (rank, start_ticks, program_name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: int64, output_unit
        use mpi_f08, only: mpi_finalize, mpi_success
        use run_timer, only: finish_run_timer
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: rank
        integer(int64), intent(in) :: start_ticks
        character(len=*), intent(in) :: program_name

        ! Local declarations
        integer :: mpi_return_code

        flush (output_unit)
        call mpi_finalize (mpi_return_code)
        if (mpi_return_code /= mpi_success) error stop "MPI_Finalize failed"
        if (rank == 0) call finish_run_timer (start_ticks, program_name)

        return
    end subroutine finish_mpi_application

    ! Print messages on rank 0 in rank order. Workers send the length before the payload.
    subroutine print_mpi_message (message, rank, process_count)

        ! Module dependencies
        use mpi_f08, only: mpi_character, mpi_comm_world, mpi_integer, mpi_recv, mpi_send, mpi_status_ignore
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  ::  rank, process_count
        character(len=*), intent(inout) :: message

        ! Local declarations
        integer :: source_rank, mpi_return_code, message_length

        message_length = len_trim (message)

        if (rank /= 0) then
            call mpi_send (message_length, 1, mpi_integer, 0, rank, &
                mpi_comm_world, mpi_return_code)
            call mpi_send (message, message_length, mpi_character, 0, rank+process_count, &
                mpi_comm_world, mpi_return_code)
        else
            print *, message(1:message_length)
            do source_rank = 1, process_count - 1
                call mpi_recv (message_length, 1, mpi_integer, source_rank, &
                    source_rank, mpi_comm_world, mpi_status_ignore, mpi_return_code)
                call mpi_recv (message, message_length, mpi_character, source_rank, &
                    source_rank+process_count, mpi_comm_world, mpi_status_ignore, mpi_return_code)
                print *, message(1:message_length)
            end do
        end if

        return
    end subroutine print_mpi_message

    ! Pause all ranks until rank 0 receives input.
    subroutine check_mpi_checkpoint (rank, checkpoint_message)
        ! Module dependencies
        use mpi_f08, only: mpi_barrier, mpi_comm_world
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: rank
        character(len=*), intent(in) :: checkpoint_message

        ! Local declarations
        integer :: mpi_return_code
        character(len=1) :: confirmation_character
        if (rank == 0) then
            print *, checkpoint_message, "hit any key to continue"
            read (*, "(a)") confirmation_character
        end if

        call mpi_barrier (mpi_comm_world, mpi_return_code)

        return
    end subroutine check_mpi_checkpoint

    ! Sum real values element by element and return the result on every rank.
    subroutine allreduce_mpi_real_sum (values, value_count)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use mpi_f08, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_sum
        use blas_lapack_interfaces, only: dcopy
        use real_vector_initialization, only: initialize_real_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in)                         :: value_count
        real(real64), intent(inout) ::   values(1:value_count)

        ! Local declarations
        integer                                     :: mpi_return_code
        real(real64)                :: summed_values(1:value_count)

        call initialize_real_vector (value_count, 0.0_real64, summed_values, 1)
        call mpi_allreduce (values, summed_values, value_count, mpi_double_precision, &
            mpi_sum, mpi_comm_world, mpi_return_code)
        call dcopy (value_count, summed_values, 1, values, 1)

        return
    end subroutine allreduce_mpi_real_sum

    ! Sum real values element by element and return the result on rank 0.
    subroutine reduce_mpi_real_sum_to_root (values, value_count)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use mpi_state, only: mpi_rank
        use mpi_f08, only: mpi_comm_world, mpi_double_precision, mpi_reduce, mpi_sum
        use blas_lapack_interfaces, only: dcopy
        use real_vector_initialization, only: initialize_real_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in)                          :: value_count
        real(real64), intent(inout) :: values(1:value_count)

        ! Local declarations
        integer                                      :: mpi_return_code
        real(real64)                 :: summed_values(1:value_count)

        if (mpi_rank == 0) then
            call initialize_real_vector (value_count, 0.0_real64, summed_values, 1)
            call mpi_reduce (values, summed_values, value_count, mpi_double_precision, mpi_sum, 0, &
                mpi_comm_world, mpi_return_code)
            call dcopy (value_count, summed_values, 1, values, 1)
        else
            call mpi_reduce (values, summed_values, value_count, mpi_double_precision, mpi_sum, 0, &
                mpi_comm_world, mpi_return_code)
        end if

        return
    end subroutine reduce_mpi_real_sum_to_root

    ! Take element-wise maxima of real values and return the result on rank 0.
    subroutine reduce_mpi_real_max_to_root (values, value_count)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use mpi_state, only: mpi_rank
        use mpi_f08, only: mpi_comm_world, mpi_double_precision, mpi_max, mpi_reduce
        use blas_lapack_interfaces, only: dcopy
        implicit none (type, external)

        ! Arguments
        integer, intent(in)                          :: value_count
        real(real64), intent(inout) :: values(1:value_count)

        ! Local declarations
        integer                                      :: mpi_return_code
        real(real64)                 :: maximum_values(1:value_count)

        if (mpi_rank == 0) then
            call mpi_reduce (values, maximum_values, value_count, mpi_double_precision, mpi_max, 0, &
                mpi_comm_world, mpi_return_code)
            call dcopy (value_count, maximum_values, 1, values, 1)
        else
            call mpi_reduce (values, maximum_values, value_count, mpi_double_precision, mpi_max, 0, &
                mpi_comm_world, mpi_return_code)
        end if

        return
    end subroutine reduce_mpi_real_max_to_root

    ! Sum integer values element by element and return the result on every rank.
    subroutine allreduce_mpi_integer_sum (values, value_count)
        ! Module dependencies
        use mpi_f08, only: mpi_allreduce, mpi_comm_world, mpi_integer, mpi_sum
        use integer_vector_copy, only: copy_integer_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in)                         :: value_count
        integer, intent(inout) ::  values(1:value_count)

        ! Local declarations
        integer                                     :: mpi_return_code
        integer                ::  summed_values(1:value_count)

        summed_values = 0

        call mpi_allreduce (values, summed_values, value_count, mpi_integer, mpi_sum, &
            mpi_comm_world, mpi_return_code)
        call copy_integer_vector (value_count, summed_values, 1, values, 1)

        return
    end subroutine allreduce_mpi_integer_sum
end module mpi_utilities
