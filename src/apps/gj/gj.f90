!***********************************************************************
!                                                                      *
program gj
    ! Module dependencies
    !                                                                      *
    !   Entry routine for the serial Landé-factor program.                 *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
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
    character :: name*24
    !-----------------------------------------------
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "GJ"
    write (stderr_unit, *) "This is the serial Lande-factor program"
    write (stderr_unit, *) "Input files:  isodata, name_desc.parquet, name_header.toml, name.mix/name.cimix, name.wfn"
    write (stderr_unit, *) "Output files: name.gj or name.cgj"
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "Default settings?"
    prompt_answer = read_yes_no()
    write (stderr_unit, *)
    if (prompt_answer) then
        nondefault_options_flag = 0
    else
        nondefault_options_flag = 1
    end if
    !
    call read_state_name(name)
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
    prompt_answer = read_yes_no()
    if (prompt_answer) then
        mixing_file_kind = 0
    else
        mixing_file_kind = 1
    end if
    !
    call configure_debug_output
    call initialize_machine_constants
    call initialize_physical_constants
    call load_angular_csf_basis(name, unused_core_subshell_count)
    call read_hyperfine_options(name)
    call read_mixing_blocks(name, mixing_file_kind)
    call initialize_log_factorial_table
    !
    call alloc(lande_g_factors, active_level_count, "GJC_DIAG", "GJ")
    call alloc(lande_g_factor_corrections, active_level_count, "DGJC_DIAG", "GJ")
    !
    call compute_gj_only(lande_g_factors, lande_g_factor_corrections)
    call write_lande_g_factors(name, mixing_file_kind, lande_g_factors, lande_g_factor_corrections)
    !
    call dalloc(lande_g_factors, "GJC_DIAG", "GJ")
    call dalloc(lande_g_factor_corrections, "DGJC_DIAG", "GJ")
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "GJ: Execution complete."
    stop
end program gj
