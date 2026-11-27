! Assemble and diagonalize the relativistic CI Hamiltonian across MPI ranks.
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
program rci_mpi_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use calculation_file_units, only: calculation_log_unit
    use calculation_input_defaults, only: nondefault_options_flag
    use rci_prerun_settings, only: prerun_stage
    use standard_io_units, only: stderr_unit
    use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_process_count, mpi_rank, processor_name, &
        processor_name_length
    use hamiltonian_correction_flags, only: include_transverse_interaction
    use calculation_state_files, only: orbital_input_path, descriptor_input_path
    use mpi_calculation_state, only: read_mpi_state_name
    use yes_no_input, only: read_yes_no
    use mpi_csf_block_input, only: open_csf_blocks_mpi
    use mpi_io_errors, only: abort_input_mpi
    use breit_radial_integral_coefficients, only: set_breit_abort_handler

    ! Procedure interfaces
    use configure_ci_debug_output_interface, only: configure_ci_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file
    use connect_ci_restart_file_interface, only: connect_ci_restart_file
    use read_ci_options_interface, only: write_ci_restart_header
    use open_mixing_output_mpi_interface, only: open_mixing_output_mpi
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use solve_hamiltonian_blocks_interface, only: solve_hamiltonian_blocks
    use tabulate_breit_integrals_1_mpi_interface, only: tabulate_breit_integrals_1_mpi
    use tabulate_breit_integrals_2_mpi_interface, only: tabulate_breit_integrals_2_mpi
    use tabulate_slater_integrals_mpi_interface, only: tabulate_slater_integrals_mpi

    ! Module dependencies
    use mpi_utilities, only: initialize_mpi_runtime_with_status
    use mpi_utilities, only: finish_mpi_application
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50
    integer, parameter :: nblk=50
    integer(int64) :: start_clock_count
    integer :: core_subshell_count, max_orbital_two_j
    logical :: yes
    logical :: is_restart
    character(len=128) :: startup_directory, state_name, scratch_directory, permanent_directory
    character(:), allocatable :: name, isotope_file_path, orbital_path, descriptor_file_path
    character(len=8) :: block_labels(maximum_input_symmetry_blocks)
    ! Initialize MPI and set up the I/O and scratch directories.

    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        & start_clock_count,       &
        startup_directory, permanent_directory, scratch_directory, "RCI_MPI")
    call set_breit_abort_handler(abort_input_mpi)

    ! Read the calculation options on rank 0 and broadcast them.

    if (mpi_rank == 0) then
        write (stderr_unit,*)
        write (stderr_unit, *) "Default settings? "
        yes = read_yes_no()
        if (yes) then
            nondefault_options_flag = 0

        else
            nondefault_options_flag = 1

        end if
    end if
    call mpi_bcast (nondefault_options_flag,1,mpi_integer,0,mpi_comm_world,mpi_error_code)

    ! Read the state name used to resolve the input and output paths.
    ! info used only on node-0

    name = ""
    isotope_file_path = ""
    orbital_path = ""
    descriptor_file_path = ""
    call read_mpi_state_name(state_name)
    if (mpi_rank == 0) then
        open(newunit=calculation_log_unit,file=trim(state_name)//".cilog",status="unknown")
        write(calculation_log_unit,"(a)") "y            ! Default settings"
        write(calculation_log_unit,"(a)") trim(state_name)
        isotope_file_path = trim(permanent_directory)//"/isodata"
        name = trim(permanent_directory)//"/"//trim(state_name)
        orbital_path = orbital_input_path(name)
        descriptor_file_path = descriptor_input_path(name)
    end if

    do

        prerun_stage = 0

        ! Open the descriptor basis and restart matrix before setting up the requested eigenpairs.

        call configure_ci_debug_output
        call initialize_machine_constants
        call initialize_physical_constants
        if (mpi_rank == 0)  call open_summary_file (name)
        call open_csf_blocks_mpi (descriptor_file_path, core_subshell_count, block_labels)
        call connect_ci_restart_file (state_name, isotope_file_path, orbital_path, block_labels, is_restart)
        ! Validate the requested eigenpairs before writing mixing and restart records.
        call open_mixing_output_mpi (name, block_labels)
        if (.not. is_restart) call write_ci_restart_header()
        if (mpi_rank == 0)  call write_input_summary
        call initialize_log_factorial_table

        ! Compute each rank’s part of the Slater integrals and assemble the shared table.

        call tabulate_slater_integrals_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)

        ! Tabulate the two Breit integral groups when that interaction is enabled.
        if (include_transverse_interaction) then
            call tabulate_breit_integrals_1_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)
            call tabulate_breit_integrals_2_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)
        end if

        ! Proceed with the CI calculation
        write (output_unit, *) "Calling solve_hamiltonian_blocks..."

        call solve_hamiltonian_blocks (core_subshell_count, max_orbital_two_j)

        if (prerun_stage /= 1) exit
        prerun_stage = 2
    end do
    ! Report elapsed time.

    call finish_mpi_application (mpi_rank, start_clock_count, "RCI_MPI")

    stop
end program rci_mpi_main
