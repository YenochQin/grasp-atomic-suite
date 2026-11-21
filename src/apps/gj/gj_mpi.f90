!***********************************************************************
!                                                                      *
program gj_mpi
    ! Module dependencies
    !                                                                      *
    !   Entry routine for the MPI Landé-factor program.                    *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
    use array_memory_management, only: alloc, dalloc
    use mpi_state, only: mpi_bcast, mpi_character, mpi_comm_world, mpi_error_code, mpi_finalize, mpi_integer, &
        mpi_process_count, mpi_rank, processor_name, processor_name_length
    use selected_level_indices,          only: active_level_count
    use compute_gj_only_mpi_interface, only: compute_gj_only_mpi

    ! Procedure interfaces
    use read_hyperfine_options_interface, only: read_hyperfine_options
    use read_mixing_blocks_interface, only: read_mixing_blocks

    ! Module dependencies
    use calculation_state_files, only: read_state_name, property_state_arguments
    use yes_no_input, only: read_yes_no
    use physical_constants, only: initialize_physical_constants
    use mpi_csf_basis_input, only: load_csf_basis_mpi

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants

    ! Procedure interfaces
    use write_lande_g_factors_interface, only: write_lande_g_factors

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table
    use mpi_utilities, only: mpix_shutdown
    use mpi_utilities, only: mpix_startup
    implicit none (type, external)
    integer :: input_status, mixing_file_kind, unused_core_subshell_count, start_clock_count
    logical :: prompt_answer, use_command_line_arguments
    real(real64), pointer :: lande_g_factors(:), lande_g_factor_corrections(:)
    character :: name*24
    !-----------------------------------------------
    !
    call mpix_startup(mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count, &
        "GJ_MPI", "This is the MPI Lande-factor program", &
        "isodata, descriptor Parquet/header, name.mix/name.cimix, name.wfn", &
        "name.gj or name.cgj", .true.)
    !
    name = " "
    mixing_file_kind = 1
    nondefault_options_flag = 0
    use_command_line_arguments = .false.
    !
    if (mpi_rank == 0) then
        use_command_line_arguments = command_argument_count() > 0
        if (use_command_line_arguments) then
            call property_state_arguments("gj_mpi", name, mixing_file_kind, input_status)
            if (input_status /= 0) nondefault_options_flag = -input_status
        end if
        if (.not.use_command_line_arguments .and. nondefault_options_flag == 0) then
            write (stderr_unit, *)
            write (stderr_unit, *) "Default settings?"
            prompt_answer = read_yes_no()
            write (stderr_unit, *)
            if (prompt_answer) then
                nondefault_options_flag = 0
            else
                nondefault_options_flag = 1
                write (stderr_unit, *) "GJ_MPI currently supports default settings only."
            end if
        end if
        !
        if (nondefault_options_flag == 0 .and. .not.use_command_line_arguments) then
            call read_state_name(name, status=input_status)
            if (input_status /= 0) then
                nondefault_options_flag = -input_status
            else
                !
                write (stderr_unit, *)
                write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
                prompt_answer = read_yes_no()
                if (prompt_answer) then
                    mixing_file_kind = 0
                else
                    mixing_file_kind = 1
                end if
            end if
        end if
    end if
    !
    call mpi_bcast(nondefault_options_flag, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    if (nondefault_options_flag /= 0) then
        call mpi_finalize(mpi_error_code)
        if (nondefault_options_flag == -2) stop 2
        stop
    end if
    !
    call mpi_bcast(name, len(name), mpi_character, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast(mixing_file_kind, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
    !
    call configure_debug_output
    call initialize_machine_constants
    call initialize_physical_constants
    call load_csf_basis_mpi(name, unused_core_subshell_count)
    call read_hyperfine_options(name)
    call read_mixing_blocks(name, mixing_file_kind)
    call initialize_log_factorial_table
    !
    call alloc(lande_g_factors, active_level_count, "lande_g_factors", "GJ_MPI")
    call alloc(lande_g_factor_corrections, active_level_count, "lande_g_factor_corrections", "GJ_MPI")
    !
    call compute_gj_only_mpi(lande_g_factors, lande_g_factor_corrections)
    if (mpi_rank == 0) call write_lande_g_factors(name, mixing_file_kind, lande_g_factors, lande_g_factor_corrections)
    !
    call dalloc(lande_g_factors, "lande_g_factors", "GJ_MPI")
    call dalloc(lande_g_factor_corrections, "lande_g_factor_corrections", "GJ_MPI")
    !
    call mpix_shutdown(mpi_rank, start_clock_count, "GJ_MPI")
    stop
end program gj_mpi
