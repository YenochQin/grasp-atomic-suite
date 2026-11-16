program cache_io
    use angdata_i
    use input_io, only: fail_input_io
    use iso_fortran_env, only: real64
#ifdef CACHE_MPI
    use mpi
    use input_io_mpi, only: fail_input_io_mpi
#endif
    implicit none (type, external)
    character(24) :: name = "cache"
    character :: argument
    logical :: available, opened
    integer :: status, rank, i, word, whichone, unit
    real(real64) :: coefficient
#ifdef CACHE_MPI
    integer :: error
#endif
    rank = 0
    call get_command_argument(1, argument)
    read(argument,*) whichone
    unit = 49 + whichone
#ifdef CACHE_MPI
    call mpi_init(error)
    call mpi_comm_rank(mpi_comm_world, rank, error)
#endif
    do i = 1, 2
        if (rank == 0) then
#ifdef CACHE_BASELINE
            call angdata(name, available, whichone)
            status = 0
#else
            call angdata(name, available, whichone, status)
#endif
            if (status /= 0) then
#ifdef CACHE_MPI
                call fail_input_io_mpi(status, unit, "ris4: angular cache")
#else
                call fail_input_io(status, unit, "ris4: angular cache")
#endif
            end if
            print *, "AVAILABLE", available
            inquire(unit=unit, opened=opened)
            if (available) then
                if (.not. opened) error stop "available cache not connected"
                read(unit) word
                if (word /= -1 .and. word /= 42) error stop "cache not rewound"
                if (word == 42) then
                    backspace(unit)
                    read(unit) word, coefficient
                    if (coefficient /= 0.125_real64) error stop "coefficient changed"
                end if
                close(unit)
#ifndef CACHE_BASELINE
            else
                if (opened) error stop "unavailable cache still connected"
#endif
            end if
        end if
#ifdef CACHE_MPI
        call mpi_bcast(available, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, error)
#endif
    end do
#ifdef CACHE_MPI
    call mpi_finalize(error)
#endif
end program cache_io
