! Calculate and report radiative multipole transitions from biorthogonal wavefunctions.
! Use orbitals transformed by rbiotransform so the two orbital sets are biorthogonal.
! Written by Per Jonsson,   Department of Computer Science
!                         Vanderbilt University, USA
!
! Modified by Gediminas Gaigalas for new spin-angular integration
! and for reducing usage of CPU memory.        NIST, October 2017
program rtransition_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use calculation_input_defaults, only: input_mixing_file_kind, nondefault_options_flag, angular_output_dump_flag
    use debug_output_settings, only: rmcdhf_debug_flags
    use yes_no_input, only: read_yes_no
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use read_state_filenames_interface, only: read_state_filenames
    use merge_transition_csf_bases_interface, only: merge_transition_csf_bases

    use read_transition_options_interface, only: read_transition_options
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use calculate_oscillator_strengths_interface, only: calculate_oscillator_strengths

    ! Module dependencies
    use run_timer, only: begin_run_timer, finish_run_timer
    implicit none (type, external)
    integer :: ntest, label_output_option
    integer(int64) :: start_clock_count
    logical :: yes
    character :: name(2)*24
    character(len=128) :: isotope_file_path

    call begin_run_timer (start_clock_count, "RTRANSITION")
    isotope_file_path = "isodata"
    ntest = 1001
    write(*,*) "This program computes transition parameters from"
    write(*,*) "transformed wave functions"
    write(*,*) "Input files:  isodata, name1_desc.parquet, name1.bw, name1.(c)bm"
    write(*,*) "              name2_desc.parquet, name2.bw, name2.(c)bm         "
    write(*,*) "              optional, name1.lsj.lbl, name2.lsj.lbl"
    write(*,*) "              name1.name2.KT (optional angular files)"
    write(*,*) "Output files: name1.name2.(c)t                       "
    write(*,*) "              optional, name1.name2.(c)t.lsj         "
    write(*,*) "              name1.name2.KT (angular files)         "
    write(*,*) "Here K is parity and rank of transition: -1,+1 etc   "

    write (output_unit, *)
    write (output_unit, *) " Default settings?"
    yes = read_yes_no()
    write (output_unit, *)
    if (yes) then
        nondefault_options_flag = 0
        angular_output_dump_flag = 1
    else
        nondefault_options_flag = 1
        write (output_unit, *) " Dump angular data to file?"
        yes = read_yes_no()
        if (yes) then
            angular_output_dump_flag = 1
        else
            angular_output_dump_flag = 0
        end if
    end if
    write (output_unit, *)
    write (output_unit, *) " Input from a CI calculation?"
    yes = read_yes_no()
    write (output_unit, *)
    if (yes) then
        input_mixing_file_kind = 0
    else
        input_mixing_file_kind = 1
    end if
    ! start
    rmcdhf_debug_flags = .false.
    ! end
    !
    ! Initialize machine-dependent constants.
    call initialize_machine_constants
    ! Initialize physical constants.
    call initialize_physical_constants
    ! Obtain the names of the initial and final state files
    call read_state_filenames (name)
    ! Load and align the initial and final CSF bases in memory.
    call merge_transition_csf_bases (name)
    !
    !
    !
    !
    !   Get the remaining information
    call read_transition_options (isotope_file_path,name)
    ! Open the transition report and record the input settings.
    label_output_option = 0
    call write_input_summary (name, input_mixing_file_kind,label_output_option)
    ! Initialize the log-factorial table for angular algebra.
    call initialize_log_factorial_table
    ! Proceed with the transition calculation
    call calculate_oscillator_strengths (name)

    call finish_run_timer (start_clock_count, "RTRANSITION")
    stop
end program rtransition_main
