program cpath_io
    use mpi_c
    use cpath_I, only: cpath
    use sys_chdir_I, only: sys_chdir
    use sys_getwd_I, only: sys_getwd
    implicit none (type, external)
    character(128) :: startdir, permdir, tmpdir
    character(128) :: expected, current
    character(3) :: rank_text
    integer :: status, length
    call startmpi(mpi_rank, mpi_process_count, processor_name, processor_name_length)
    call get_environment_variable("MPI_TMP", expected)
    write(rank_text,"(I3.3)") mpi_rank
    call cpath(startdir, permdir, tmpdir)
    if (trim(tmpdir) /= trim(expected)) error stop "CPATH temporary root"
    call sys_getwd(current, length)
    if (trim(current) /= trim(expected)//"/"//rank_text) error stop "CPATH rank directory"
    if (mpi_rank == 0) then
        if (permdir /= startdir) error stop "CPATH permanent directory"
    end if
    call sys_chdir(startdir, len_trim(startdir), status)
    if (status /= 0) error stop "CPATH restore working directory"
    call MPI_Finalize(status)
    print *, "PASS CPATH input environment"
end program cpath_io
