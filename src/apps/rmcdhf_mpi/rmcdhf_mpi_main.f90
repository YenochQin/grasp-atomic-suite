! Optimize orbitals with the MPI MCDHF solver using rank-local MCP files.
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
program rmcdhf_mpi_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64
    use calculation_file_units, only: calculation_log_unit
    use calculation_input_defaults, only: nondefault_options_flag
    use closed_core_subshells, only: closed_core_subshell_count
    use mpi_state, only: mpi_process_count, mpi_rank, processor_name, processor_name_length
    use calculation_state_files, only: check_orbital_output
    use calculation_state_files, only: state_files, state_name_length, files_for_state
    use mpi_calculation_state, only: select_mpi_state

    ! Procedure interfaces
    use configure_debug_output_mpi_interface, only: configure_debug_output_mpi

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file
    use connect_mcp_files_interface, only: connect_mcp_files
    use load_csf_basis_data_mpi_interface, only: load_csf_basis_data_mpi
    use read_scf_options_mpi_interface, only: read_scf_options_mpi
    use write_input_summary_interface, only: write_input_summary
    use open_mixing_output_interface, only: open_mixing_output

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use run_scf_iterations_mpi_interface, only: run_scf_iterations_mpi

    ! Module dependencies
    use mpi_utilities, only: initialize_mpi_runtime_with_status
    use mpi_utilities, only: finish_mpi_application
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50  ! Maximum number of blocks
    integer :: core_subshell_count, lenperm
    integer(int64) :: start_clock_count
    type(state_files) :: files
    character(state_name_length) :: state
    logical :: extended_optimal_level, yes
    character :: block_labels(maximum_input_symmetry_blocks)*8
    character (len = 128) :: startup_directory,serial_io_directory,rank_scratch_root
    character(:), allocatable :: csf_descriptor_path, isotope_file_path, initial_radial_file_path
    character(:), allocatable :: mixing_output_path, radial_output_path

    ! Initialize MPI and set up the I/O and scratch directories.
    startup_directory = "  "
    csf_descriptor_path = "  "
    serial_io_directory = "  "
    isotope_file_path = "  "
    rank_scratch_root = "  "
    initial_radial_file_path = "  "
    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        & start_clock_count,       &
        startup_directory, serial_io_directory, rank_scratch_root, "RMCDHF_MPI")
    lenperm = len_trim (serial_io_directory)

    ! Record whether the run uses default SCF options.

    call select_mpi_state(state, "Default settings? (y/n)", yes, serial_io_directory)
    files = files_for_state(state)
    if (mpi_rank == 0) call check_orbital_output(trim(serial_io_directory)//"/"//files%orbitals, writing_initial=.false.)
    nondefault_options_flag = merge(0, 1, yes)
    if (mpi_rank == 0) then
        open(newunit=calculation_log_unit,file=trim(startup_directory)//"/"//files%scf_log,status="UNKNOWN")
        if (state /= "") then
            write(calculation_log_unit,"(A)") trim(state)
            write(calculation_log_unit,"(A)") merge("y","n",yes)//"            ! Default settings"
        end if
    end if

    ! Load the MCP metadata and the descriptor basis, then check that they agree.
    call configure_debug_output_mpi (trim(serial_io_directory)//"/"//files%scf_debug)
    call initialize_machine_constants
    call initialize_physical_constants

    if (mpi_rank == 0) call open_summary_file (trim(serial_io_directory)//"/"//files%scf_summary)

    call connect_mcp_files (closed_core_subshell_count, maximum_input_symmetry_blocks, block_labels, files%mcp)
    if(mpi_rank == 0) csf_descriptor_path = trim(serial_io_directory)//"/"//files%descriptors
    call load_csf_basis_data_mpi (csf_descriptor_path, core_subshell_count, block_labels)
    if (closed_core_subshell_count /= core_subshell_count) error stop "rmcdhf_mpi_main: closed-core subshell count mismatch"

    ! Read the SCF controls and prepare the radial orbitals.

    if(mpi_rank == 0) isotope_file_path = serial_io_directory(1:lenperm) // "/isodata"
    if(mpi_rank == 0) initial_radial_file_path = trim(serial_io_directory)//"/"//files%initial_orbitals
    call read_scf_options_mpi (extended_optimal_level, block_labels, isotope_file_path, initial_radial_file_path )
    mixing_output_path = "  "
    if(mpi_rank == 0) mixing_output_path = trim(serial_io_directory)//"/"//files%mixing
    if (mpi_rank == 0) then
        call write_input_summary
        if (extended_optimal_level) call open_mixing_output (mixing_output_path)
    end if

    call initialize_log_factorial_table

    ! Iterate the SCF solution and save the optimized orbitals.

    radial_output_path = "  "
    if(mpi_rank == 0) radial_output_path = trim(serial_io_directory)//"/"//files%orbitals

    call run_scf_iterations_mpi (extended_optimal_level, radial_output_path)

    ! Report elapsed time.

    call finish_mpi_application (mpi_rank, start_clock_count, "RMCDHF_MPI")

end program rmcdhf_mpi_main
