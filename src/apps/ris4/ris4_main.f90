! Calculate the electronic mass-shift and field-shift factors used in isotope shifts.
! Written by Per Jonsson. Last revision: 17 Jan 1996
!    Modify  by Gediminas Gaigalas                        26 Oct 2009
! Modified by J. Ekman. 25 Nov 2013
!
! Translated by Gediminas Gaigalas 11/18/19
program ris4_main
    !  Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
    use run_timer, only: begin_run_timer, finish_run_timer
    use debug_output_settings,         only: matrix_element_screening_cutoff
    use calculation_state_files, only: read_state_name
    use yes_no_input, only: read_yes_no

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file

    ! Module dependencies
    use csf_basis_loader, only: load_angular_csf_basis

    ! Procedure interfaces
    use read_mixing_blocks_interface, only: read_mixing_blocks
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use calculate_isotope_shifts_interface, only: calculate_isotope_shifts
    use prepare_isotope_calculation_interface, only: prepare_isotope_calculation
    implicit none (type, external)
    logical   :: yes
    character :: name*24
    integer :: mixing_file_kind, unused_core_subshell_count
    integer(int64) :: start_clock_count

    call begin_run_timer (start_clock_count, "ris4_main", log_unit=stderr_unit)
    ! Matrix elements smaller than CUTOFF are not accumulated
    matrix_element_screening_cutoff = 1.0e-10_real64
    write (stderr_unit, *)
    write (stderr_unit, *) "Default settings?"
    yes = read_yes_no ()
    write (stderr_unit, *)
    if (yes) then
        nondefault_options_flag = 0
    else
        nondefault_options_flag = 1
    end if

    call read_state_name(name)
    print *

    write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
    yes = read_yes_no ()
    if (yes) then
        mixing_file_kind = 0
    else
        mixing_file_kind = 1
    end if
    write (stderr_unit, *)

    call configure_debug_output
    ! Initialize machine-dependent constants.
    call initialize_machine_constants
    ! Initialize physical constants.
    call initialize_physical_constants
    ! Open the  .i  file
    call open_summary_file(name,mixing_file_kind)
    ! Load the configuration state functions from the descriptor pair.
    call load_angular_csf_basis(name,unused_core_subshell_count)
    ! Get the remaining information
    call prepare_isotope_calculation(name)
    ! Get the eigenvectors

    call read_mixing_blocks(name,mixing_file_kind)
    ! Record the input settings in the calculation report.
    !
    ! Initialize the log-factorial table for angular algebra.
    call initialize_log_factorial_table
    ! Proceed with the RIS calculation
    call calculate_isotope_shifts(name)

    call finish_run_timer (start_clock_count, "ris4_main", log_unit=stderr_unit)
    stop
end  program ris4_main
