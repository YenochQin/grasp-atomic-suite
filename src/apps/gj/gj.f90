!***********************************************************************
!                                                                      *
program gj
    ! Module dependencies
    !                                                                      *
    !   Entry routine for the serial Landé-factor program.                 *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use calculation_input_defaults, only: inpci, ndef, ndump
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use array_memory_management, only: alloc, dalloc
    use selected_level_indices,          only: nvec
    use compute_gj_only_interface, only: compute_gj_only

    ! Procedure interfaces
    use read_hyperfine_options_interface, only: read_hyperfine_options
    use read_mixing_blocks_interface, only: read_mixing_blocks

    ! Module dependencies
    use calculation_state_files, only: read_state_name
    use yes_no_input, only: getyn
    use physical_constants, only: setcon
    use csf_basis_loader, only: setcsla

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: setmc

    ! Procedure interfaces
    use write_lande_g_factors_interface, only: write_lande_g_factors

    ! Module dependencies
    use log_factorial_table, only: factt
    implicit none (type, external)
    integer :: nci, unused_core_subshell_count
    logical :: yes
    real(real64), pointer :: gjc_diag(:), dgjc_diag(:)
    character :: name*24
    !-----------------------------------------------
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "GJ"
    write (stderr_unit, *) "This is the serial Lande-factor program"
    write (stderr_unit, *) "Input files:  isodata, name_descriptors.parquet, name_header.toml, name.mix/name.cimix, name.wfn"
    write (stderr_unit, *) "Output files: name.gj or name.cgj"
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "Default settings?"
    yes = getyn()
    write (stderr_unit, *)
    if (yes) then
        ndef = 0
    else
        ndef = 1
    end if
    !
    call read_state_name(name)
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
    yes = getyn()
    if (yes) then
        nci = 0
    else
        nci = 1
    end if
    !
    call configure_debug_output
    call setmc
    call setcon
    call setcsla(name, unused_core_subshell_count)
    call read_hyperfine_options(name)
    call read_mixing_blocks(name, nci)
    call factt
    !
    call alloc(gjc_diag, nvec, "GJC_DIAG", "GJ")
    call alloc(dgjc_diag, nvec, "DGJC_DIAG", "GJ")
    !
    call compute_gj_only(gjc_diag, dgjc_diag)
    call write_lande_g_factors(name, nci, gjc_diag, dgjc_diag)
    !
    call dalloc(gjc_diag, "GJC_DIAG", "GJ")
    call dalloc(dgjc_diag, "DGJC_DIAG", "GJ")
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "GJ: Execution complete."
    stop
end program gj
