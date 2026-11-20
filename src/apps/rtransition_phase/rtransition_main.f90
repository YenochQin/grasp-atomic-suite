!***********************************************************************
!                                                                      *
program rtransition_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !     This program calculates the transition parameters for a          *
    !     transition between an initial and a final state                  *
    !     The program assumes that the radial orbitals of the initial      *
    !     and final state have been transformed by the BIOTRA program      *
    !     as to become biorthonormal, in which case the normal Racah       *
    !     algebra can be used.                                             *
    !                                                                      *
    !     Written by Per Jonsson,   Department of Computer Science         *
    !                             Vanderbilt University, USA               *
    !                                                                      *
    !     Modified by Gediminas Gaigalas for new spin-angular integration  *
    !     and for reducing usage of CPU memory.        NIST, October 2017  *
    !                                                                      *
    !***********************************************************************
    use calculation_input_defaults, only: inpci, ndef, ndump
    use debug_output_settings, only: ldbpr, cutoff
    use yes_no_input, only: getyn
    use machine_constants, only: setmc
    use physical_constants, only: setcon

    ! Procedure interfaces
    use read_state_filenames_interface, only: read_state_filenames
    use merge_transition_csf_bases_interface, only: merge_transition_csf_bases

    use read_transition_options_interface, only: read_transition_options
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: factt

    ! Procedure interfaces
    use calculate_oscillator_strengths_interface, only: calculate_oscillator_strengths

    ! Module dependencies
    use run_timer_start, only: starttime
    use run_timer_finish, only: stoptime
    implicit none (type, external)
    integer :: ntest, start_clock_count, ilbl
    logical :: yes
    character :: name(2)*24
    character(len=128) :: isotope_file_path
    !-----------------------------------------------
    !
    isotope_file_path = "isodata"
    ntest = 1001
    call starttime (start_clock_count, "RTRANSITION")
    write(*,*)
    write(*,*) "RTRANSITION"
    write(*,*) "This program computes transition parameters from"
    write(*,*) "transformed wave functions"
    write(*,*) "Input files:  isodata, name1_descriptors.parquet, name1.bw, name1.(c)bm"
    write(*,*) "              name2_descriptors.parquet, name2.bw, name2.(c)bm         "
    write(*,*) "              optional, name1.lsj.lbl, name2.lsj.lbl"
    write(*,*) "              name1.name2.KT (optional angular files)"
    write(*,*) "Output files: name1.name2.(c)t                       "
    write(*,*) "              optional, name1.name2.(c)t.lsj         "
    write(*,*) "              name1.name2.KT (angular files)         "
    write(*,*) "Here K is parity and rank of transition: -1,+1 etc   "




    write (output_unit, *)
    write (output_unit, *) " Default settings?"
    yes = getyn()
    write (output_unit, *)
    if (yes) then
        ndef = 0
        ndump = 1
    else
        ndef = 1
        write (output_unit, *) " Dump angular data to file?"
        yes = getyn()
        if (yes) then
            ndump = 1
        else
            ndump = 0
        end if
    end if
    write (output_unit, *)
    write (output_unit, *) " Input from a CI calculation?"
    yes = getyn()
    write (output_unit, *)
    if (yes) then
        inpci = 0
    else
        inpci = 1
    end if
    !Rasa -- start
    ldbpr = .false.
    !      WRITE (6, *) ' Generate debug output?'
    !      YES = GETYN()
    !      WRITE (6, *)
    !      IF (YES) THEN
    !         LDBPR(18) = .TRUE.
    !         WRITE (6, *) ' Enter the cutoff'
    !         READ (5, *) CUTOFF
    !      ENDIF
    !Rasa -- end
    !
    !   Perform machine- and installation-dependent setup
    !
    call setmc
    !
    !   Set up the physical constants
    !
    call setcon
    !
    !   Obtain the names of the initial and final state files
    !
    call read_state_filenames (name)
    !
    !   Open, check, load data from, and close, the initial and final state
    !   .csl  files. These files are then merged to one file.
    !
    call merge_transition_csf_bases (name)
    !
    !   Open, check, load data from, and close, the merged .csl  file
    !
    ! The descriptor bases were merged directly in memory by MRGCSL.
    !
    !   Read mixing coefficients
    !
    !      CALL READMIX(NAME,INPCI)
    !
    !   Test mixing coefficients
    !
    !      IF (NTEST.GT.1000) CALL TESTMIX
    !
    !   Get the remaining information
    !
    call read_transition_options (isotope_file_path,name)
    !
    !   Open and append a summary of the inputs to the  .sum  file
    !
    ilbl = 0
    call write_input_summary (name, inpci,ilbl)
    !
    !   Set up the table of logarithms of factorials
    !
    call factt
    !
    !   Proceed with the transition calculation
    !
    call calculate_oscillator_strengths (name)
    !
    !   Print completion message
    !
    call stoptime (start_clock_count, "RTRANSITION")
    !
    stop
end program rtransition_main
