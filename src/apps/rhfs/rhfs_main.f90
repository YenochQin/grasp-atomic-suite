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
program rhfs_main
    ! Module dependencies
    !                                                                      *
    !   Entry routine for HFS92. Controls the entire computation.          *
    !                                                                      *
    !   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
    !               [HFS92]: CHKPLT, GETHFD, HFS, SETDBG, SETSUM,          *
    !                        STRSUM.                                       *
    !               [NJGRAF]: FACTT.                                       *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
    !   Modified by Gediminas Gaigalas for new spin-angular integration.   *
    !                                         Last revision:    Nov 2017   *
    !                                                                      *
    !***********************************************************************
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
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
    !-----------------------------------------------
    !

    write (stderr_unit, *)
    write (stderr_unit, *) "RHFS"
    write (stderr_unit, *) "This is the hyperfine structure program"
    write (stderr_unit, *) "Input files:  isodata, name_desc.parquet, name_header.toml, name.mix/name.cimix, name.wfn"
    write (stderr_unit, *) "Output files: name.(c)h, name.(c)hoffd"

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
    !
    !   Check compatibility of plant substitutions
    !
    !GG      CALL CHKPLT
    !
    !   Determine if there is to be any debug printout; this will be
    !   made on the  .dbg  file
    !
    call configure_debug_output
    !
    !   Perform machine- and installation-dependent setup
    !
    call initialize_machine_constants
    !
    !   Set up the physical constants
    !
    call initialize_physical_constants
    !
    !   Open the  .sum  file
    !
    call open_summary_file (name, mixing_file_kind)
    !
    !   Open, check, load data from, and close, the  .csl  file
    !
    call load_angular_csf_basis (name, unused_core_subshell_count)
    !
    !   Get the remaining information
    !
    call read_hyperfine_options (name)
    !
    !   Get the eigenvectors
    !
    !GG      WRITE (ISTDE, *) 'Block format?'
    !GG      YES = GETYN()
    !GG      WRITE (ISTDE, *)
    !GG      IF (YES) THEN
    call read_mixing_blocks (name, mixing_file_kind)
    !GG      ELSE
    !GG         IF (NCI == 0) THEN
    !GG            CALL GETMIXC (NAME)
    !GG         ELSE
    !GG            CALL GETMIXA (NAME)
    !GG         ENDIF
    !GG      ENDIF
    !
    !   Append a summary of the inputs to the  .sum  file
    !
    call write_input_summary
    !
    !   Set up the table of logarithms of factorials
    !
    call initialize_log_factorial_table
    !
    !   Proceed with the HFS calculation
    !
    call calculate_hyperfine_matrix_elements
    !
    !   Print completion message
    !
    write (stderr_unit, *)
    write (stderr_unit, *) "RHFS: Execution complete."
    !
    stop
end program rhfs_main
