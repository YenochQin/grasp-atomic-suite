! Calculate hyperfine constants and Landé factors for selected relativistic states.
!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**           **   **  *******   *****    *****    *****             ***
!**           **   **  **       **   **  **   **  **   **            ***
!**           **   **  **       **       **   **      **             ***
!**           *******  ****      *****    *****      **              ***
!**           **   **  **            **      **     **               ***
!**           **   **  **       **   **     **     **                ***
!**           **   **  **        *****    **      *******            ***
!**                                                                  ***
!**            Relativistic Hyperfine Structure Program              ***
!**                         GRASP92 Version                          ***
!**                         Dynamic Storage                          ***
!**                                                                  ***
!**   ============================================================   ***
!**   Copyright (c) 1995 by P Jonsson, F A Parpia, and C F Fischer   ***
!**   ============================================================   ***
!**   All rights reserved.  No part of this software or its accom-   ***
!**   panying documentation may be reproduced, stored in a retrie-   ***
!**   val system,  or transmitted,  in any form or  by any  means,   ***
!**   electronic, mechanical,  photocopying,  recording, or other-   ***
!**   wise, without the prior written permission of the authors.     ***
!**                                                                  ***
!**                           Disclaimer                             ***
!**                           ==========                             ***
!**   The  authors make  no warranties,  express or implied,  that   ***
!**   this software or its  accompanying documentation are free of   ***
!**   error or that  they will meet your requirements for any par-   ***
!**   ticular application.  Neither this software nor its accompa-   ***
!**   nying documentation should be relied  upon for solving prob-   ***
!**   lems if an incorrect solution could result in injury or loss   ***
!**   of, or damage to, property. If you use this software or  its   ***
!**   accompanying documentation,  you do so entirely  at your own   ***
!**   risk;  the authors disclaim all liability for direct or con-   ***
!**   sequential damage.                                             ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
! Written by Farid A. Parpia. Last revision: 24 Dec 1992
! Modified by Gediminas Gaigalas for new spin-angular integration.
! Last revision: Nov 2017
program hyperfine_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
    use run_timer, only: begin_run_timer, finish_run_timer
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
    use read_hyperfine_options_interface, only: read_hyperfine_options
    use read_mixing_blocks_interface, only: read_mixing_blocks
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use calculate_hyperfine_matrix_elements_interface, only: calculate_hyperfine_matrix_elements
    implicit none (type, external)
    integer   :: mixing_file_kind, unused_core_subshell_count
    logical   :: prompt_answer
    character :: name*24
    integer(int64) :: start_clock_count

    call begin_run_timer (start_clock_count, "HYPERFINE", log_unit=stderr_unit)
    write (stderr_unit, *) "This is the hyperfine structure program"
    write (stderr_unit, *) "Input files:  isodata, name_desc.parquet, name_header.toml, name.mix/name.cimix, name.wfn"
    write (stderr_unit, *) "Output files: name.(ci)hyperfine, name.(ci)hyperfine.off-diagonal"

    write (stderr_unit, *)
    write (stderr_unit, *) "Default settings?"
    prompt_answer = read_yes_no()
    write (stderr_unit, *)
    if (prompt_answer) then
        nondefault_options_flag = 0
    else
        nondefault_options_flag = 1
    end if

    call read_state_name(name)
    write (stderr_unit, *)
    write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
    prompt_answer = read_yes_no()
    if (prompt_answer) then
        mixing_file_kind = 0
    else
        mixing_file_kind = 1
    end if
    call configure_debug_output("hyperfine.dbg")
    ! Initialize machine-dependent constants.
    call initialize_machine_constants
    ! Initialize physical constants.
    call initialize_physical_constants
    ! Open the  .sum  file
    call open_summary_file (name, mixing_file_kind)
    ! Load the configuration state functions from the descriptor pair.
    call load_angular_csf_basis (name, unused_core_subshell_count)
    ! Get the remaining information
    call read_hyperfine_options (name)
    !    Get the eigenvectors
    call read_mixing_blocks (name, mixing_file_kind)
    ! Record the input settings in the calculation report.
    call write_input_summary
    ! Initialize the log-factorial table for angular algebra.
    call initialize_log_factorial_table
    ! Proceed with the HFS calculation
    call calculate_hyperfine_matrix_elements

    call finish_run_timer (start_clock_count, "HYPERFINE", log_unit=stderr_unit)
    stop
end program hyperfine_main
