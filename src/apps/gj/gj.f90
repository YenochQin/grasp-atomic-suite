! Calculate electronic Landé factors and anomalous magnetic moment corrections for selected levels.
program gj
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
    use run_timer, only: begin_run_timer, finish_run_timer
    use array_memory_management, only: alloc, dalloc
    use selected_level_indices,          only: active_level_count
    use compute_gj_only_interface, only: compute_gj_only

    ! Procedure interfaces
    use read_hyperfine_options_interface, only: read_hyperfine_options
    use read_mixing_blocks_interface, only: read_mixing_blocks

    ! Module dependencies
    use calculation_state_files, only: read_state_name
    use yes_no_input, only: read_yes_no
    use physical_constants, only: initialize_physical_constants
    use csf_basis_loader, only: load_angular_csf_basis

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants

    ! Procedure interfaces
    use write_lande_g_factors_interface, only: write_lande_g_factors

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table
    implicit none (type, external)
    integer :: mixing_file_kind, unused_core_subshell_count
    logical :: prompt_answer
    real(real64), pointer :: lande_g_factors(:), lande_g_factor_corrections(:)
    character :: state_name*24
    integer(int64) :: start_clock_count

    call begin_run_timer (start_clock_count, "GJ", log_unit=stderr_unit)
    write (stderr_unit, *) "This is the serial Lande-factor program"
    write (stderr_unit, *) "Input files:  isodata, name_desc.parquet, name_header.toml, name.mix/name.cimix, name.wfn"
    write (stderr_unit, *) "Output files: name.gj or name.cigj"
    write (stderr_unit, *)
    write (stderr_unit, *) "Default settings?"
    prompt_answer = read_yes_no()
    write (stderr_unit, *)
    if (prompt_answer) then
        nondefault_options_flag = 0
    else
        nondefault_options_flag = 1
    end if
    call read_state_name(state_name)
    write (stderr_unit, *)
    write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
    prompt_answer = read_yes_no()
    if (prompt_answer) then
        mixing_file_kind = 0
    else
        mixing_file_kind = 1
    end if
    call configure_debug_output("gj.dbg")
    call initialize_machine_constants
    call initialize_physical_constants
    call load_angular_csf_basis(state_name, unused_core_subshell_count)
    call read_hyperfine_options(state_name)
    call read_mixing_blocks(state_name, mixing_file_kind)
    call initialize_log_factorial_table
    call alloc(lande_g_factors, active_level_count, "lande_g_factors", "GJ")
    call alloc(lande_g_factor_corrections, active_level_count, "lande_g_factor_corrections", "GJ")
    call compute_gj_only(lande_g_factors, lande_g_factor_corrections)
    call write_lande_g_factors(state_name, mixing_file_kind, lande_g_factors, lande_g_factor_corrections)
    call dalloc(lande_g_factors, "lande_g_factors", "GJ")
    call dalloc(lande_g_factor_corrections, "lande_g_factor_corrections", "GJ")
    call finish_run_timer (start_clock_count, "GJ", log_unit=stderr_unit)
    stop
end program gj
