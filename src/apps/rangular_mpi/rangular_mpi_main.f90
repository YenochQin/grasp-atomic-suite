! Generate rank-local MCP channels for the MPI MCDHF solvers.
! Program for generating the energy expression for H(DC). This
! program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,
! and C. F. Fischer, 1990).
!
! GRASP92 Version
! F. A. Parpia, C. F. Fischer, and I. P. Grant
! Written by Farid A. Parpia. Last revision: 11 Dec 1992
! MPI version by Xinghong He. Last revision: 29 Jun 1998
! Updated by Charlotte F. Fischer
!
! Modified by Gediminas Gaigalas for new spin-angular integration.
program rangular_mpi_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    use calculation_file_units, only: angular_log_unit, summary_output_unit
    use debug_output_settings, only: angular_debug_unit
    use mpi_csf_collectives, only: allocate_descriptor_arrays
    use iso_fortran_env, only: int64
    use rangular_mpi_input, only: mpi_angular_input_options, parse_mpi_angular_options, &
        open_mpi_angular_profile, load_mpi_angular_header, load_mpi_angular_block, close_mpi_angular_input
    use rangular_profile, only: stage_clock, record_stage, close_stage_profile
    use mpi_state, only: mpi_process_count, mpi_rank, processor_name, processor_name_length
    use calculation_input_defaults, only: nondefault_options_flag
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations
    use csf_coupling_data, only: csf_angular_codes, csf_subshell_quantum_codes
    use mcp_generation_options, only: max_coulomb_multipole_rank
    use calculation_state_files, only: state_files, state_name_length, files_for_state, resolve_descriptor_input
    use mpi_calculation_state, only: select_mpi_state

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file

    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use open_temporary_mcp_channels_interface, only: open_temporary_mcp_channels

    ! Module dependencies
    use mpi_utilities, only: initialize_mpi_runtime_with_status
    use mpi_utilities, only: finish_mpi_application

    ! Procedure interfaces
    use configure_debug_output_mpi_interface, only: configure_debug_output_mpi
    use open_mcp_outputs_mpi_interface, only: open_mcp_outputs_mpi
    use compute_mcp_angular_coefficients_mpi_interface, only: compute_mcp_angular_coefficients_mpi
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50   !-----------------------------------------------
    integer :: core_subshell_count, symmetry_block_index
    integer(int64) :: start_clock_count
    logical ::  debug, is_restart, yes
    character(len=128) :: startup_directory, permanent_directory, scratch_directory
    character(len=8):: block_labels(maximum_input_symmetry_blocks)
    character(len=3) :: block_label_text
    type(mpi_angular_input_options) :: input_options
    type(state_files) :: files
    character(state_name_length) :: state
    integer(int64) :: program_start_count, stage_start_count
    character(len=16) :: input_format_name
    call parse_mpi_angular_options(input_options)
    startup_directory = " "
    permanent_directory = " "
    scratch_directory = " "
    ! Initialize MPI and set up the I/O and scratch directories.
    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        & start_clock_count, &
        startup_directory, permanent_directory, scratch_directory, "RANGULAR_MPI")
    write (block_label_text, "(I3.3)") mpi_rank
    call select_mpi_state(state, "Use zero-first partition from descriptor header? (y/n)", yes, permanent_directory)
    files = files_for_state(state)
    if (input_options%input_filename == "") then
        input_options%input_filename = files%descriptors
        if (mpi_rank == 0) input_options%input_filename = resolve_descriptor_input(files%descriptors, permanent_directory)
    end if
    if (mpi_rank == 0) then
        open(newunit=angular_log_unit,file=trim(startup_directory)//"/"//files%angular_log,status="UNKNOWN")
        if (state /= "") write(angular_log_unit,"(A)") trim(state)
    end if
    call open_mpi_angular_profile(input_options)
    input_format_name = "descriptors"
    program_start_count = stage_clock()
    ! Log whether to take the zero-first partition from the descriptor header.
    nondefault_options_flag = merge(1, 0, yes)
    if (mpi_rank == 0) write(angular_log_unit,"(A)") &
        merge("y", "n", yes)//"            ! Use zero-first partition from descriptor header"

    ! Read the descriptor header and initialize the output channels for the angular coefficients.

    call configure_debug_output_mpi (debug, trim(permanent_directory)//"/"//files%angular_debug)
    call initialize_machine_constants
    stage_start_count = stage_clock()
    call load_mpi_angular_header(input_options,permanent_directory,core_subshell_count,maximum_input_symmetry_blocks,block_labels)
    call record_stage("input_header_and_broadcast",0,stage_start_count,storage=input_format_name)
    stage_start_count = stage_clock()
    is_restart = .false.

    call open_mcp_outputs_mpi (core_subshell_count, block_labels, files%mcp, input_options%descriptors)
    if (nondefault_options_flag/=0 .and. mpi_rank==0) call write_input_summary
    call initialize_log_factorial_table
    call record_stage("mcp_setup_including_interaction_input",0,stage_start_count)

    ! Generate and sort the angular contributions one symmetry block at a time.

    do symmetry_block_index = 1, symmetry_block_count
        csf_count = csf_counts_by_block(symmetry_block_index)                        ! This ncf goes to common
        if (mpi_rank == 0) then
            write (output_unit, *)
            write (output_unit, *) "Block ", symmetry_block_index, ",  ncf = ", csf_count
        end if
        ! Load the current descriptor block before computing its angular coefficients.
        stage_start_count = stage_clock()
        call allocate_descriptor_arrays(csf_count,csf_subshell_occupations, &
            csf_subshell_quantum_codes,csf_angular_codes)
        call record_stage("csf_allocation",symmetry_block_index,stage_start_count,int(csf_count,int64))
        stage_start_count = stage_clock()
        call load_mpi_angular_block(input_options,core_subshell_count,symmetry_block_index)
        call record_stage("csf_load_and_broadcast",symmetry_block_index,stage_start_count,int(csf_count,int64), &
            & storage=input_format_name)
        ! Open the temporary MCP channels for this symmetry block.
        call open_temporary_mcp_channels (symmetry_block_index, max_coulomb_multipole_rank, files%temporary//block_label_text)
        ! Generate this symmetry block’s MCP data.
        call compute_mcp_angular_coefficients_mpi (symmetry_block_index, is_restart, files%mcp)
    end do
    call close_mpi_angular_input(input_options)
    if (summary_output_unit /= -1) close(summary_output_unit)                                  ! Summary file
    if (mpi_rank == 0) close(angular_log_unit)                                 ! rangular.log
    if (debug .and. angular_debug_unit /= -1) close(angular_debug_unit)                       ! Debug file
    call record_stage("program_before_finalize_including_interaction_input",0,program_start_count)
    call close_stage_profile()
    ! Report elapsed time.
    call finish_mpi_application (mpi_rank, start_clock_count, "RANGULAR_MPI")
    stop
end program rangular_mpi_main
