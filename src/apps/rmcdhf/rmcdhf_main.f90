! Optimize radial orbitals and CSF mixing coefficients with a serial MCDHF calculation.
! Relativistic Self-Consistent-Field Program
!
! This program is a derivative of GRASP2 (F. A. Parpia, I. P.
! Grant, and C. F. Fischer, 1990)
!
! GRASP92
! F. A. Parpia, C. F. Fischer, and I. P. Grant
! Modified by C.F. FIscher for block input
! Written by Farid A. Parpia. Last revision: 31 Dec 1992
! Block version by Xinghong He. Last revision: 17 Aug 1998
! Modified by G. Gaigalas. 05 Feb 2017
program rmcdhf_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64
    use calculation_file_units, only: calculation_log_unit
    use calculation_input_defaults, only: nondefault_options_flag
    use closed_core_subshells, only: closed_core_subshell_count
    use serial_mpi_context, only: mpi_process_count, mpi_rank
    use calculation_state_files, only: check_orbital_output
    use calculation_state_files, only: state_files, state_name_length, files_for_state, select_state
    use run_timer, only: begin_run_timer, finish_run_timer

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file
    use connect_mcp_files_interface, only: connect_mcp_files
    use load_csf_basis_data_interface, only: load_csf_basis_data
    use read_scf_options_interface, only: read_scf_options
    use write_input_summary_interface, only: write_input_summary
    use open_mixing_output_interface, only: open_mixing_output

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use run_scf_iterations_interface, only: run_scf_iterations

    ! Module dependencies
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50  ! Maximum number of blocks
    integer :: core_subshell_count
    integer(int64) :: start_clock_count
    type(state_files) :: files
    character(state_name_length) :: state
    logical :: extended_optimal_level, yes
    character :: block_labels(maximum_input_symmetry_blocks)*8
    call begin_run_timer (start_clock_count, "RMCDHF")
    ! Use a single-process context for the serial kernels.
    mpi_rank = 0
    mpi_process_count = 1

    write(*,*) "This program determines the radial orbitals   "
    write(*,*) "and the expansion coefficients of the CSFs         "
    write(*,*) "in a self-onsistent field proceedure               "
    call select_state(state, "Default settings? (y/n)", yes)
    files = files_for_state(state)
    call check_orbital_output(files%orbitals, writing_initial=.false.)
    write(*,*) "Input files: isodata, ", files%descriptors, ", ", files%header, ", ", files%initial_orbitals
    write(*,*) "Output files: ", files%orbitals, ", ", files%mixing, ", ", files%scf_summary, ", ", files%scf_log


    open(newunit=calculation_log_unit,file=files%scf_log,status="UNKNOWN")
    if (state /= "") write(calculation_log_unit,"(A)") trim(state)

    ! Record whether the run uses default SCF options.

    if (mpi_rank == 0) then
        if (yes) then
            nondefault_options_flag = 0
            write(calculation_log_unit,"(A)") "y            ! Default settings"
        else
            nondefault_options_flag = 1
            if (state /= "") write(calculation_log_unit,"(A)") "n            ! Default settings"
        end if
    end if

    ! Load the MCP metadata and the descriptor basis, then check that they agree.

    call configure_debug_output (files%scf_debug)
    call initialize_machine_constants
    call initialize_physical_constants

    call open_summary_file (files%scf_summary)

    call connect_mcp_files (closed_core_subshell_count, maximum_input_symmetry_blocks, block_labels, files%mcp)
    call load_csf_basis_data (files%descriptors, core_subshell_count, block_labels)
    if (closed_core_subshell_count /= core_subshell_count) error stop "RMCDHF: closed-core subshell count mismatch"

    ! Read the SCF controls and prepare the radial orbitals.

    call read_scf_options (extended_optimal_level, block_labels, "isodata", files%initial_orbitals, "npol.inp")

    if (mpi_rank == 0) then
        call write_input_summary
        if (extended_optimal_level) call open_mixing_output (files%mixing)
    end if

    call initialize_log_factorial_table

    ! Iterate the SCF solution and save the optimized orbitals.

    call run_scf_iterations (extended_optimal_level, files%orbitals)
    close (calculation_log_unit)

    ! Report elapsed time.

    call finish_run_timer (start_clock_count, "RMCDHF")

    stop
end program rmcdhf_main
