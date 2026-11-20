module mpi_utilities
    ! Module dependencies
    use mpi_runtime, only: startmpi, mpix_bytes, stopmpi
    implicit none (type, external)

    ! Module interface
    private
    public :: startmpi
    public :: mpix_startup
    public :: startmpi2
    public :: mpix_shutdown
    public :: stopmpi
    public :: stopmpi2
    public :: mpix_printmsg
    public :: mpix_chkpt
    public :: mpix_bytes
    public :: gdsummpi
    public :: gdrsummpi_root
    public :: gdmaxmpi_root
    public :: gisummpi
contains

    ! Start MPI in the caller's working directory and optionally redirect worker output.
    subroutine mpix_startup (rank, process_count, processor_name, processor_name_length, start_ticks, &
        program_name, description, input_file_description, &
        output_file_description, quiet_workers)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: error_unit

        implicit none (type, external)

        ! Arguments
        integer :: rank, process_count, processor_name_length, start_ticks
        character(len=*), intent(out) :: processor_name
        character(len=*), intent(in) :: program_name, description
        character(len=*), intent(in) :: input_file_description, output_file_description
        logical, intent(in) :: quiet_workers

        ! Local declarations
        integer :: ticks_per_second, clock_count_limit, output_open_status
        character(len=256) :: output_message

        call startmpi (rank, process_count, processor_name, processor_name_length)
        call system_clock (start_ticks, ticks_per_second, clock_count_limit)

        if (rank == 0) then
            write (output_unit, *)
            write (output_unit, *) trim(program_name)
            write (output_unit, *) trim(description)
            write (output_unit, *) "Input files:  " // trim(input_file_description)
            write (output_unit, *) "Output files: " // trim(output_file_description)
            write (output_unit, *) "Running on ", process_count, " MPI ranks"
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

        return
    end subroutine mpix_startup

    ! Start MPI, report the participating ranks and set up rank-local working directories.
    subroutine startmpi2 (rank, process_count, processor_name, processor_name_length, start_ticks, &
        startup_directory, serial_io_directory, rank_scratch_root, program_name)

        ! Module dependencies
        use mpi_work_paths, only: cpath

        implicit none (type, external)

        ! Arguments
        integer, intent(out)  :: rank,  process_count, processor_name_length
        integer, intent(out)  :: start_ticks
        character(len=*), intent(out) :: processor_name
        character(len=*), intent(in) :: program_name
        character(len=*), intent(out) :: startup_directory, serial_io_directory, rank_scratch_root

        ! Local declarations
        character :: rank_label*3
        character :: date_text*8, time_text*10, zone_text*5
        integer ::  date_time_values(8)
        integer :: startup_directory_length, serial_io_directory_length, rank_scratch_root_length
        integer :: ticks_per_second, clock_count_limit
        character(len=80) :: message

        call startmpi (rank, process_count, processor_name, processor_name_length)
        write (rank_label, "(I3.3)") rank
        if (rank == 0) then
            print *, "===================================================="
            print *, "       ", program_name, ": Execution Begins ..."
            print *, "===================================================="
            print *,        "Participating nodes:"
        end if
        message = "  Host: " // processor_name(1:processor_name_length) // "    ID: " // rank_label
        call mpix_printmsg (message, rank, process_count)

        call date_and_time (date_text, time_text, zone_text, date_time_values)
        if (rank == 0) then
            print *
            print *, "Date and Time:"
        end if
        message = "  " // processor_name(1:processor_name_length) // ": " // &
            "  Date: " // date_text // &
            "  Time: " // time_text // &
            "  Zone: " // zone_text
        call mpix_printmsg (message, rank, process_count)

        ! Keep the serial I/O directory for rank 0 and select a scratch root for each rank.
        call cpath (startup_directory, serial_io_directory, rank_scratch_root)

        startup_directory_length = len_trim (startup_directory)
        serial_io_directory_length = len_trim (serial_io_directory)
        rank_scratch_root_length = len_trim (rank_scratch_root)

        if (rank == 0) then
            print *
            print *, "Start Dir:"
        end if
        message = "  " // processor_name(1:processor_name_length) // ": " // startup_directory(1:startup_directory_length)
        call mpix_printmsg (message, rank, process_count)

        if (rank == 0) then
            print *
            print *, "Serial I/O Dir (node-0 only):"
            print *, "  " // processor_name(1:processor_name_length) // ": " // serial_io_directory(1:serial_io_directory_length)
        end if

        if (rank == 0) then
            print *
            print *, "Work Dir (Parallel I/O):"
        end if
        message = "  " // processor_name(1:processor_name_length) // ": " // rank_scratch_root(1:rank_scratch_root_length)
        call mpix_printmsg (message, rank, process_count)

        call system_clock (start_ticks, ticks_per_second, clock_count_limit)
        return
    end subroutine startmpi2

    ! Print elapsed time on rank 0 and finalize MPI.
    subroutine mpix_shutdown (rank, start_ticks, program_name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use mpi, only: mpi_finalize
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: rank, start_ticks
        character(len=*), intent(in) :: program_name

        ! Local declarations
        integer :: elapsed_ticks, ticks_per_second, clock_count_limit, elapsed_seconds, mpi_return_code
        character(len=8) :: date_text
        character(len=10) :: time_text
        character(len=5) :: zone_text
        integer :: date_time_values(8)

        call system_clock (elapsed_ticks, ticks_per_second, clock_count_limit)
        elapsed_ticks = elapsed_ticks - start_ticks
        elapsed_seconds = elapsed_ticks / ticks_per_second
        call date_and_time (date_text, time_text, zone_text, date_time_values)

        if (rank == 0) then
            write (output_unit, *)
            write (output_unit, *) trim(program_name) // ": Execution complete."
            write (output_unit, *) "Wall time: ", elapsed_seconds, " seconds"
            write (output_unit, *) "Finish date: ", date_text, " time: ", time_text, &
                " zone: ", zone_text
        end if

        call mpi_finalize (mpi_return_code)

        return
    end subroutine mpix_shutdown

    ! Collect elapsed time and completion details on rank 0, then synchronize and stop.
    subroutine stopmpi2 (rank, process_count, processor_name, processor_name_length, start_ticks, &
        program_name)

        ! Module dependencies
        use mpi, only: mpi_barrier, mpi_comm_world

        implicit none (type, external)

        ! Arguments
        integer, intent(in)  ::  rank,process_count,processor_name_length
        integer              ::  start_ticks
        character*(*), intent(in) ::  processor_name, program_name
        integer ::   elapsed_ticks, ticks_per_second, clock_count_limit, elapsed_seconds
        character :: date_text*8, time_text*10, zone_text*5
        integer ::  date_time_values(8)

        character :: elapsed_seconds_text*8, message*80

        integer :: mpi_return_code

        if (rank == 0) then
            print *, "===================================================="
            print *, "       ", program_name, ": Execution Finished ..."
            print *, "===================================================="
            print *,        "Wall time:"
        end if

        call system_clock (elapsed_ticks, ticks_per_second, clock_count_limit)
        elapsed_ticks = elapsed_ticks - start_ticks
        elapsed_seconds = elapsed_ticks / ticks_per_second
        write (elapsed_seconds_text, "(i8)") elapsed_seconds
        message = elapsed_seconds_text // " seconds on " // processor_name(1:processor_name_length)
        call mpix_printmsg (message, rank, process_count)

        if (rank == 0) then
            print *
            print *, "Finish Date and Time:"
        end if

        call date_and_time (date_text, time_text, zone_text, date_time_values)

        message = "  " // processor_name(1:processor_name_length) // ": " // &
            "  Date: " // date_text // &
            "  Time: " // time_text // &
            "  Zone: " // zone_text
        call mpix_printmsg (message, rank, process_count)

        if (rank == 0) print *

        call mpi_barrier (mpi_comm_world,mpi_return_code)
        call stopmpi (program_name // ": Execution complete.", rank)
        return
    end subroutine stopmpi2

    ! Print messages on rank 0 in rank order. Workers send the length before the payload.
    subroutine mpix_printmsg (message, rank, process_count)

        ! Module dependencies
        use mpi, only: mpi_character, mpi_comm_world, mpi_integer, mpi_recv, mpi_send, mpi_status_size
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  ::  rank, process_count
        character(len=*), intent(inout) :: message

        ! Local declarations
        integer :: source_rank, receive_status(mpi_status_size), mpi_return_code, message_length

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
                    source_rank, mpi_comm_world, receive_status, mpi_return_code)
                call mpi_recv (message, message_length, mpi_character, source_rank, &
                    source_rank+process_count, mpi_comm_world, receive_status, mpi_return_code)
                print *, message(1:message_length)
            end do
        end if

        return
    end subroutine mpix_printmsg

    ! Pause all ranks until rank 0 receives input.
    subroutine mpix_chkpt (rank, checkpoint_message)
        ! Module dependencies
        use mpi, only: mpi_barrier, mpi_comm_world
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
    end subroutine mpix_chkpt

    ! Sum real values element by element and return the result on every rank.
    subroutine gdsummpi (values, value_count)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use mpi, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_sum
        use blas_lapack_interfaces, only: dcopy
        use real_vector_initialization, only: dinit
        implicit none (type, external)

        ! Arguments
        integer, intent(in)                         :: value_count
        real(real64), intent(inout) ::   values(1:value_count)

        ! Local declarations
        integer                                     :: mpi_return_code
        real(real64)                :: summed_values(1:value_count)

        call dinit (value_count, 0.0_real64, summed_values, 1)
        call mpi_allreduce (values, summed_values, value_count, mpi_double_precision, &
            mpi_sum, mpi_comm_world, mpi_return_code)
        call dcopy (value_count, summed_values, 1, values, 1)

        return
    end subroutine gdsummpi

    ! Sum real values element by element and return the result on rank 0.
    subroutine gdrsummpi_root (values, value_count)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use mpi_state, only: mpi_rank
        use mpi, only: mpi_comm_world, mpi_double_precision, mpi_reduce, mpi_sum
        use blas_lapack_interfaces, only: dcopy
        use real_vector_initialization, only: dinit
        implicit none (type, external)

        ! Arguments
        integer, intent(in)                          :: value_count
        real(real64), intent(inout) :: values(1:value_count)

        ! Local declarations
        integer                                      :: mpi_return_code
        real(real64)                 :: summed_values(1:value_count)

        if (mpi_rank == 0) then
            call dinit (value_count, 0.0_real64, summed_values, 1)
            call mpi_reduce (values, summed_values, value_count, mpi_double_precision, mpi_sum, 0, &
                mpi_comm_world, mpi_return_code)
            call dcopy (value_count, summed_values, 1, values, 1)
        else
            call mpi_reduce (values, summed_values, value_count, mpi_double_precision, mpi_sum, 0, &
                mpi_comm_world, mpi_return_code)
        end if

        return
    end subroutine gdrsummpi_root

    ! Take element-wise maxima of real values and return the result on rank 0.
    subroutine gdmaxmpi_root (values, value_count)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use mpi_state, only: mpi_rank
        use mpi, only: mpi_comm_world, mpi_double_precision, mpi_max, mpi_reduce
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
    end subroutine gdmaxmpi_root

    ! Sum integer values element by element and return the result on every rank.
    subroutine gisummpi (values, value_count)
        ! Module dependencies
        use mpi, only: mpi_allreduce, mpi_comm_world, mpi_integer, mpi_sum
        use integer_vector_copy, only: icopy
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
        call icopy (value_count, summed_values, 1, values, 1)

        return
    end subroutine gisummpi
end module mpi_utilities
