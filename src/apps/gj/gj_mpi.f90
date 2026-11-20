!***********************************************************************
!                                                                      *
program gj_mpi
    ! Module dependencies
    !                                                                      *
    !   Entry routine for the MPI Landé-factor program.                    *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use calculation_input_defaults, only: inpci, ndef, ndump
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use array_memory_management, only: alloc, dalloc
    use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
        & mpi_process_count, mpi_rank, processor_name, processor_name_length
    use selected_level_indices,          only: nvec
    use compute_gj_only_mpi_interface, only: compute_gj_only_mpi

    ! Procedure interfaces
    use read_hyperfine_options_interface, only: read_hyperfine_options
    use read_mixing_blocks_interface, only: read_mixing_blocks

    ! Module dependencies
    use calculation_state_files, only: read_state_name, property_state_arguments
    use yes_no_input, only: getyn
    use physical_constants, only: setcon
    use mpi_csf_basis_input, only: setcsla => load_csf_basis_mpi

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: setmc

    ! Procedure interfaces
    use write_lande_g_factors_interface, only: write_lande_g_factors

    ! Module dependencies
    use log_factorial_table, only: factt
    use mpi_utilities, only: mpix_shutdown
    use mpi_utilities, only: mpix_startup
    implicit none (type, external)
    integer :: input_status, nci, unused_core_subshell_count, start_clock_count
    logical :: yes, use_args
    real(real64), pointer :: gjc_diag(:), dgjc_diag(:)
    character :: name*24
    !-----------------------------------------------
    !
    call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count, &
        "GJ_MPI", "This is the MPI Lande-factor program", &
        "isodata, descriptor Parquet/header, name.mix/name.cimix, name.wfn", &
        "name.gj or name.cgj", .true.)
    !
    name = " "
    nci = 1
    ndef = 0
    use_args = .false.
    !
    if (mpi_rank == 0) then
        use_args = command_argument_count() > 0
        if (use_args) then
            call property_state_arguments("gj_mpi", name, nci, input_status)
            if (input_status /= 0) ndef = -input_status
        end if
        if (.not.use_args .and. ndef == 0) then
            write (stderr_unit, *)
            write (stderr_unit, *) "Default settings?"
            yes = getyn()
            write (stderr_unit, *)
            if (yes) then
                ndef = 0
            else
                ndef = 1
                write (stderr_unit, *) "GJ_MPI currently supports default settings only."
            end if
        end if
        !
        if (ndef == 0 .and. .not.use_args) then
            call read_state_name(name, status=input_status)
            if (input_status /= 0) then
                ndef = -input_status
            else
                !
                write (stderr_unit, *)
                write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
                yes = getyn()
                if (yes) then
                    nci = 0
                else
                    nci = 1
                end if
            end if
        end if
    end if
    !
    call mpi_bcast(ndef, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    if (ndef /= 0) then
        call mpi_finalize(mpi_error_code)
        if (ndef == -2) stop 2
        stop
    end if
    !
    call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast(nci, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    !
    call configure_debug_output
    call setmc
    call setcon
    call setcsla(name, unused_core_subshell_count)
    call read_hyperfine_options(name)
    call read_mixing_blocks(name, nci)
    call factt
    !
    call alloc(gjc_diag, nvec, "GJC_DIAG", "GJ_MPI")
    call alloc(dgjc_diag, nvec, "DGJC_DIAG", "GJ_MPI")
    !
    call compute_gj_only_mpi(gjc_diag, dgjc_diag)
    if (mpi_rank == 0) call write_lande_g_factors(name, nci, gjc_diag, dgjc_diag)
    !
    call dalloc(gjc_diag, "GJC_DIAG", "GJ_MPI")
    call dalloc(dgjc_diag, "DGJC_DIAG", "GJ_MPI")
    !
    call mpix_shutdown(mpi_rank, start_clock_count, "GJ_MPI")
    stop
end program gj_mpi
