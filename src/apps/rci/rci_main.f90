!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**             ******    *****   ****   *****    *****              ***
!**             **   **  **   **   **   **   **  **   **             ***
!**             **   **  **        **   **   **       **             ***
!**             ******   **        **    *****       **              ***
!**             **  **   **        **      **       **               ***
!**             **   **  **   **   **     **      **                 ***
!**             **   **   *****   ****   **      *******             ***
!**                                                                  ***
!**          Relativistic Configuration-Interaction Program          ***
!**                                                                  ***
!**   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
!**   Grant, and C. F. Fischer, 1990).                               ***
!**                                                                  ***
!**                            GRASP92                               ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
program rci_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    use calculation_file_units, only: calculation_log_unit
    !                                                                      *
    !   Entry routine for RCI92. Controls the entire computation.          *
    !                                                                      *
    !   Call(s) to: [LIB92]: SETMC, SETCON.                                *
    !               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
    !                        SETRES, SETSUM, STRSUM.                       *
    !               [NJGRAF]: FACTT.                                       *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
    !   Updated by Xinghong He                Last revision: 23 Jun 1998   *
    !   Modified by Gediminas Gaigalas for new spin-angular integration.   *
    !                                                                      *
    !***********************************************************************
    use calculation_input_defaults, only: nondefault_options_flag
    use rci_prerun_settings, only: prerun_stage
    use standard_io_units, only: stderr_unit
    use hamiltonian_correction_flags, only: include_transverse_interaction
    use calculation_state_files, only: read_state_name, orbital_input_path, descriptor_input_path
    use yes_no_input, only: read_yes_no

    ! Procedure interfaces
    use configure_ci_debug_output_interface, only: configure_ci_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file
    use load_csf_basis_data_interface, only: load_csf_basis_data
    use connect_ci_restart_file_interface, only: connect_ci_restart_file
    use open_mixing_output_interface, only: open_mixing_output
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use solve_hamiltonian_blocks_interface, only: solve_hamiltonian_blocks
    use tabulate_breit_integrals_1_interface, only: tabulate_breit_integrals_1
    use tabulate_breit_integrals_2_interface, only: tabulate_breit_integrals_2
    use tabulate_slater_integrals_interface, only: tabulate_slater_integrals

    ! Module dependencies
    use run_timer_start, only: begin_run_timer
    use run_timer_finish, only: finish_run_timer
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50
    integer, parameter :: nblk=50
    integer :: start_clock_count
    integer :: mpi_rank, mpi_process_count, core_subshell_count, ndum, max_orbital_two_j
    logical :: yes
    character(len=128) :: name, isotope_file_path
    character(len=8) :: block_labels(maximum_input_symmetry_blocks)
    !-----------------------------------------------
    !

    prerun_stage = 0
    mpi_rank = 0
    mpi_process_count = 1

    write(*,*)
    write(*,*) "RCI"
    write(*,*) "This is the configuration interaction program "
    write(*,*) "Input file:  isodata, name_desc.parquet, name_header.toml, name.wfn"
    write(*,*) "Outputfiles: name.cimix, name.csum, name.clog "
    write(*,*) "             rci.res (can be used for restart)"
    write(*,*)


    !
    ! Start timing
    !
    call begin_run_timer (start_clock_count, "RCI")
    !
    ! Get NDEF
    !
    !      WRITE (ISTDE, *) 'RCI2: Execution begins ...'
    !      WRITE (ISTDE, *)
    write (stderr_unit, "(A)") "Default settings? "
    yes = read_yes_no()
    if (yes) then
        nondefault_options_flag = 0
    else
        nondefault_options_flag = 1
    end if
    !
    ! Get name of the state (used in files like <name>.c, <name>.s)
    !
    call read_state_name(name)

    ! Now the name of the state is known, open the log file

    open(newunit=calculation_log_unit, file=trim(name)//".clog",status="unknown")
    write(calculation_log_unit,"(a)") "y            ! Default settings"
    write(calculation_log_unit,"(a)") trim(name)

    !         ...Form the full name of the files used on node-0

    isotope_file_path = "isodata"
    !      WRITE (6, *) 'isofile = ', ISOFILE(1:LEN_TRIM(ISOFILE))
    !      WRITE (6, *) 'name = ', NAME(1:LEN_TRIM(NAME))

    do

        !
        ! In SETDBG of this version all control logicals are set to
        ! false thus no debug output will be made
        !
        !      WRITE (6, *) 'Calling SETDBG...'
        call configure_ci_debug_output
        !
        ! Perform machine- and installation-dependent setup
        !
        !      WRITE (6, *) 'Calling SETMC...'
        call initialize_machine_constants
        !
        ! Set up the physical constants
        !
        !      WRITE (6, *) 'Calling SETCON...'
        call initialize_physical_constants
        !
        ! Open summary file
        !
        !      WRITE (6, *) 'Calling SETSUM...'
        call open_summary_file (name)

        !      WRITE (6, *) 'Calling setcsl...'
        call load_csf_basis_data (descriptor_input_path(name), core_subshell_count, maximum_input_symmetry_blocks, block_labels)
        !
        ! Set up the  .res  file; determine if this is a restart.
        !
        !      WRITE (6, *) 'Calling SETRES...'
        call connect_ci_restart_file (isotope_file_path, orbital_input_path(name), block_labels)
        !
        !   Open the  .mix  file; determine the eigenpairs required
        !
        !      WRITE (6, *) 'Calling SETMIX...'
        call open_mixing_output (name, block_labels)
        !
        !   Append a summary of the inputs to the  .sum  file
        !
        write (output_unit, *) "Calling STRSUM..."
        call write_input_summary
        !
        !   Set up the table of logarithms of factorials
        !
        write (output_unit, *) "Calling FACTT..."
        call initialize_log_factorial_table
        !
        !   Calculate all the needed Rk integrals
        !
        write (output_unit, *) "Calling GENINTRK..."
        call tabulate_slater_integrals (mpi_rank, mpi_process_count, ndum, max_orbital_two_j)
        !
        !   If transverse interaction comput Breit integrals of type 1 and 2
        !
        if (include_transverse_interaction) then
            print *, "Calling GENINTBREIT1..."
            call tabulate_breit_integrals_1 (mpi_rank, mpi_process_count, ndum, max_orbital_two_j)
            print *, "Calling GENINTBREIT2..."
            call tabulate_breit_integrals_2 (mpi_rank, mpi_process_count, ndum, max_orbital_two_j)
        end if
        !
        !   Proceed with the CI calculation
        !
        write (output_unit, *) "Calling MATRIX..."

        call solve_hamiltonian_blocks (core_subshell_count, max_orbital_two_j)

        if (prerun_stage /= 1) exit
        prerun_stage = 2
    end do

    if (mpi_rank == 0) then
        write (output_unit, *)
        write (output_unit, *)
        write (output_unit, *) "Finish time, Statistics"
        write (output_unit, *)
    end if
    close(calculation_log_unit)
    call finish_run_timer (start_clock_count, "RCI")
    !      CALL SYSTEM_CLOCK (NCOUNT2, NCOUNT_RATE, NCOUNT_MAX)
    !      NCOUNT2 = NCOUNT2 - NCOUNT1
    !      NSECONDS = NCOUNT2/NCOUNT_RATE
    !      WRITE (STR, '(I8)') NSECONDS
    !      MSG = STR//' seconds '
    !      WRITE (6, *) MSG
    !
    !      CALL DATE_AND_TIME (CHDATE, CHTIME, CHZONE, NYMDUHMSM)
    !
    !      MSG = ' Date: '//CHDATE//' Time: '//CHTIME//' Zone: '//CHZONE
    !      WRITE (6, *) MSG
    !
    !   Print completion message
    !
    !      WRITE (6, *) 'RCI2: Execution complete.'
    !
    stop
end program rci_main
