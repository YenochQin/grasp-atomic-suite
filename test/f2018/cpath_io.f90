program cpath_io
    use, intrinsic :: iso_fortran_env, only: int64
    use mpi_state, only: mpi_process_count, mpi_rank, processor_name, processor_name_length
    use mpi_utilities, only: initialize_mpi_runtime_with_status, finish_mpi_application
    use mpi_change_directory, only: change_mpi_directory
    use mpi_current_directory, only: get_mpi_working_directory
    implicit none (type, external)
    character(128) :: startdir, permdir, tmpdir
    character(128) :: expected, current
    character(3) :: rank_text
    integer :: status, length
    integer(int64) :: start_ticks

    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        start_ticks, startdir, permdir, tmpdir, "CPATH_IO")
    if (command_argument_count() == 0) then
        call get_environment_variable("MPI_TMP", expected)
    else
        call get_command_argument(mpi_rank + 1, expected)
    end if
    write(rank_text,"(I3.3)") mpi_rank
    if (trim(tmpdir) /= trim(expected)) error stop "CPATH temporary root"
    call get_mpi_working_directory(current, length)
    if (trim(current) /= trim(expected)//"/"//rank_text) error stop "CPATH rank directory"
    if (mpi_rank == 0) then
        if (permdir /= startdir) error stop "CPATH permanent directory"
    end if
    call change_mpi_directory(startdir, len_trim(startdir), status)
    if (status /= 0) error stop "CPATH restore working directory"
    if (command_argument_count() > 0 .and. mpi_rank == mpi_process_count - 1) then
        call execute_command_line("sleep 2", exitstat=status)
        if (status /= 0) error stop "CPATH worker delay"
    end if
    call finish_mpi_application (mpi_rank, start_ticks, "CPATH_IO")
end program cpath_io
