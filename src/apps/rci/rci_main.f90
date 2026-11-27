! Solve the serial relativistic CI problem for a fixed set of radial orbitals.
! Relativistic Configuration-Interaction Program
!
! This program is a derivative of GRASP2 (F. A. Parpia, I. P.
! Grant, and C. F. Fischer, 1990).
!
! GRASP92
! F. A. Parpia, C. F. Fischer, and I. P. Grant
! Written by Farid A. Parpia. Last revision: 15 Oct 1992
! Updated by Xinghong He. Last revision: 23 Jun 1998
! Modified by Gediminas Gaigalas for new spin-angular integration.
program rci_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use calculation_file_units, only: calculation_log_unit
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
    use read_ci_options_interface, only: write_ci_restart_header
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
    use run_timer, only: begin_run_timer, finish_run_timer
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50
    integer, parameter :: nblk=50
    integer(int64) :: start_clock_count
    integer :: mpi_rank, mpi_process_count, core_subshell_count, integral_count, max_orbital_two_j
    logical :: yes
    logical :: is_restart
    character(len=128) :: name, isotope_file_path
    character(len=8) :: block_labels(maximum_input_symmetry_blocks)

    call begin_run_timer (start_clock_count, "RCI")
    prerun_stage = 0
    mpi_rank = 0
    mpi_process_count = 1

    write(*,*) "This is the configuration interaction program "
    write(*,*) "Input file:  isodata, name_desc.parquet, name_header.toml, name.wfn"
    write(*,*) "Outputfiles: name.cimix, name.cisum, name.cilog "
    write(*,*) "             name.rci.restart (can be used for restart)"
    write(*,*)

    ! Read the default settings and restart choices.
    write (stderr_unit, "(A)") "Default settings? "
    yes = read_yes_no()
    if (yes) then
        nondefault_options_flag = 0
    else
        nondefault_options_flag = 1
    end if
    ! Read the state name used to resolve the input and output paths.
    call read_state_name(name)

    ! Now the name of the state is known, open the log file

    open(newunit=calculation_log_unit, file=trim(name)//".cilog",status="unknown")
    write(calculation_log_unit,"(a)") "y            ! Default settings"
    write(calculation_log_unit,"(a)") trim(name)

    ! ...Form the full name of the files used on node-0

    isotope_file_path = "isodata"

    do

        call configure_ci_debug_output
        ! Initialize machine-dependent constants.
        call initialize_machine_constants
        ! Initialize physical constants.
        call initialize_physical_constants
        ! Open summary file
        call open_summary_file (name)

        call load_csf_basis_data (descriptor_input_path(name), core_subshell_count, maximum_input_symmetry_blocks, block_labels)
        ! Set up the state-specific restart file; determine if this is a restart.
        call connect_ci_restart_file (name, isotope_file_path, orbital_input_path(name), block_labels, is_restart)
        ! Validate the requested eigenpairs before writing mixing and restart records.
        call open_mixing_output (name, block_labels)
        if (.not. is_restart) call write_ci_restart_header()
        ! Record the input settings in the calculation report.
        write (output_unit, *) "Calling write_input_summary..."
        call write_input_summary
        ! Initialize the log-factorial table for angular algebra.
        write (output_unit, *) "Calling initialize_log_factorial_table..."
        call initialize_log_factorial_table
        ! Calculate all the needed Rk integrals
        write (output_unit, *) "Calling tabulate_slater_integrals..."
        call tabulate_slater_integrals (mpi_rank, mpi_process_count, integral_count, max_orbital_two_j)
        ! Tabulate the two Breit integral groups when that interaction is enabled.
        if (include_transverse_interaction) then
            print *, "Calling tabulate_breit_integrals_1..."
            call tabulate_breit_integrals_1 (mpi_rank, mpi_process_count, integral_count, max_orbital_two_j)
            print *, "Calling tabulate_breit_integrals_2..."
            call tabulate_breit_integrals_2 (mpi_rank, mpi_process_count, integral_count, max_orbital_two_j)
        end if
        ! Proceed with the CI calculation
        write (output_unit, *) "Calling solve_hamiltonian_blocks..."

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
    !    NCOUNT2 = NCOUNT2 - NCOUNT1
    !    NSECONDS = NCOUNT2/NCOUNT_RATE
    !    MSG = STR//' seconds '
    !
    !    MSG = ' Date: '//CHDATE//' Time: '//CHTIME//' Zone: '//CHZONE
    stop
end program rci_main
