program cpath_io
    use mpi_state, only: mpi_finalize, mpi_process_count, mpi_rank, processor_name, processor_name_length, initialize_mpi_runtime
    use mpi_work_paths, only: prepare_mpi_work_directories
    use mpi_change_directory, only: change_mpi_directory
    use mpi_current_directory, only: get_mpi_working_directory
    implicit none (type, external)
    character(128) :: startdir, permdir, tmpdir
    character(128) :: expected, current
    character(3) :: rank_text
    integer :: status, length
    call initialize_mpi_runtime(mpi_rank, mpi_process_count, processor_name, processor_name_length)
    call get_environment_variable("MPI_TMP", expected)
    write(rank_text,"(I3.3)") mpi_rank
    call prepare_mpi_work_directories(startdir, permdir, tmpdir)
    if (trim(tmpdir) /= trim(expected)) error stop "CPATH temporary root"
    call get_mpi_working_directory(current, length)
    if (trim(current) /= trim(expected)//"/"//rank_text) error stop "CPATH rank directory"
    if (mpi_rank == 0) then
        if (permdir /= startdir) error stop "CPATH permanent directory"
    end if
    call change_mpi_directory(startdir, len_trim(startdir), status)
    if (status /= 0) error stop "CPATH restore working directory"
    call MPI_Finalize(status)
    print *, "PASS CPATH input environment"
end program cpath_io
