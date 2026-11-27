program cache_io
    use open_isotope_shift_angular_cache_interface, only: open_isotope_shift_angular_cache
#ifndef CACHE_BASELINE
    use csf_descriptor_reader, only: loaded_csf_basis_id
    use debug_output_settings, only: matrix_element_screening_cutoff
    use open_isotope_shift_angular_cache_interface, only: finish_angular_cache
#endif
    use input_error_handling, only: fail_input_io
    use iso_fortran_env, only: real64
#ifdef CACHE_MPI
    use mpi_f08, only: mpi_bcast, mpi_comm_rank, mpi_comm_world, mpi_finalize, mpi_init, mpi_logical
    use mpi_io_errors, only: fail_input_io_mpi
#endif
    implicit none (type, external)
    character(24) :: name = "cache"
    character :: argument
    logical :: available, opened
    integer :: status, rank, i, word, whichone, cache_unit
    real(real64) :: coefficient
#ifdef CACHE_MPI
    integer :: error
#endif
    character(16) :: scenario
    rank = 0
    call get_command_argument(1, argument)
    read(argument,*) whichone
#ifdef CACHE_MPI
    call mpi_init(error)
    call mpi_comm_rank(mpi_comm_world, rank, error)
#endif
#ifndef CACHE_BASELINE
    call get_command_argument(2, scenario)
    loaded_csf_basis_id = repeat("a",64)
    matrix_element_screening_cutoff = 1e-10_real64
    if (rank == 0 .and. (scenario == "stamp" .or. scenario == "cutoff")) then
        open(newunit=cache_unit,file=merge("cache.IOB","cache.ITB",whichone==1), &
            status="old",form="unformatted",position="append")
        call finish_angular_cache(cache_unit, whichone)
        close(cache_unit)
    end if
    if (scenario == "cutoff") matrix_element_screening_cutoff = 1e-8_real64
#endif
    do i = 1, 2
        if (rank == 0) then
#ifdef CACHE_BASELINE
            call open_isotope_shift_angular_cache(name, available, whichone)
            ! The external baseline reader owns the connection.
            if (whichone == 1) then
                inquire(file="cache.IOB", number=cache_unit)
            else
                inquire(file="cache.ITB", number=cache_unit)
            end if
            status = 0
#else
            call open_isotope_shift_angular_cache(name, available, whichone, status, cache_unit)
#endif
            if (status /= 0) then
#ifdef CACHE_MPI
                call fail_input_io_mpi(status, cache_unit, "ris4: angular cache")
#else
                call fail_input_io(status, cache_unit, "ris4: angular cache")
#endif
            end if
            print *, "AVAILABLE", available
            opened = .false.
            if (cache_unit /= -1) inquire(unit=cache_unit, opened=opened)
            if (available) then
                if (.not. opened) error stop "available cache not connected"
                read(cache_unit) word
                if (word /= -1 .and. word /= 42) error stop "cache not rewound"
                if (word == 42) then
                    backspace(cache_unit)
                    read(cache_unit) word, coefficient
                    if (coefficient /= 0.125_real64) error stop "coefficient changed"
                end if
                close(cache_unit)
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
