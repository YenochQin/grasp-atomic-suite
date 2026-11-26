! Generate sparse Hamiltonian structure and angular coefficients for serial MCDHF runs.

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
! Modified by Gediminas Gaigalas: Feb 2017
!                             1) for new spin-angular integration,
!                             2) for sorting in the memory.
program rangular_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use calculation_file_units, only: angular_log_unit, summary_output_unit
    use debug_output_settings, only: angular_debug_unit
    use grasp_parameters,   only:  max_orbitals
    use rangular_input, only: angular_input_options, parse_angular_options, load_angular_header, &
        load_angular_block, close_angular_input
    use rangular_profile, only: stage_clock, record_stage, close_stage_profile
    use calculation_input_defaults, only: nondefault_options_flag
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals
    use csf_coupling_data, only: csf_angular_codes, csf_subshell_quantum_codes, csf_angular_code_count
    use calculation_state_files, only: state_files, state_name_length, files_for_state, select_state
    use run_timer, only: begin_run_timer, finish_run_timer

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file

    use open_mcp_outputs_interface, only: open_mcp_outputs
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use log_factorial_table, only: initialize_log_factorial_table

    ! Procedure interfaces
    use open_temporary_mcp_channels_interface, only: open_temporary_mcp_channels

    ! Procedure interfaces
    use compute_mcp_angular_coefficients_interface, only: compute_mcp_angular_coefficients
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50   !-----------------------------------------------
    integer(int64) :: start_clock_count
    integer :: core_subshell_count, symmetry_block_index
    logical ::  debug, is_restart, yes
    character(len=8):: block_labels(maximum_input_symmetry_blocks)
    integer :: mpi_rank,mpi_process_count
    type(angular_input_options) :: input_options
    type(state_files) :: files
    character(state_name_length) :: state
    integer(int64) :: program_start_count, stage_start_count
    character(len=16) :: input_format_name

    call parse_angular_options(input_options)
    call begin_run_timer (start_clock_count, "RANGULAR")
    call select_state(state, "Full interaction? (y/n)", yes)
    files = files_for_state(state)
    if (input_options%input_filename == "") input_options%input_filename = files%descriptors
    input_format_name = "descriptors"
    program_start_count = stage_clock()
    write(*,*) "This program performs angular integration "
    write(*,*) "Input file: ", input_options%input_filename
    write(*,*) "Output files: ", files%mcp, ".30, ", files%mcp, ".31, ..."
    write(*,*) files%angular_log
    write(*,*)
    open(newunit=angular_log_unit,file=files%angular_log,status="UNKNOWN")
    if (state /= "") write(angular_log_unit,"(A)") trim(state)
    mpi_rank = 0
    mpi_process_count = 1
    ! Record the choice of full or restricted CSF interactions.
    if (mpi_rank == 0) then
        if (yes) then
            nondefault_options_flag = 0
            write(angular_log_unit,"(A)") "y            ! Full interaction"
        else
            nondefault_options_flag = 1
            write(angular_log_unit,"(A)") "n            ! Full interaction"
        end if
    end if
    ! Read the descriptor header and initialize the output channels for the angular coefficients.
    ! delete chkplt
    call configure_debug_output (debug, files%angular_debug)
    call initialize_machine_constants
    stage_start_count = stage_clock()
    call load_angular_header(input_options, core_subshell_count, maximum_input_symmetry_blocks, block_labels)
    call record_stage("input_header", 0, stage_start_count, int(csf_count,int64), storage=input_format_name)
    is_restart = .false.
    stage_start_count = stage_clock()
    call open_mcp_outputs (mpi_rank, mpi_process_count, core_subshell_count, block_labels, files%mcp)
    if (nondefault_options_flag/=0 .and. mpi_rank==0) call write_input_summary
    call initialize_log_factorial_table
    call record_stage("setup_including_interaction_input", 0, stage_start_count)
    ! Generate and sort the angular contributions one symmetry block at a time.
    do symmetry_block_index = 1, symmetry_block_count
        csf_count = csf_counts_by_block(symmetry_block_index)                        ! This ncf goes to common
        if (mpi_rank == 0) then
            write (output_unit, *)
            write (output_unit, *) "Block ", symmetry_block_index, ",  ncf = ", csf_count
        end if
        ! Load the current descriptor block before computing its angular coefficients.
        stage_start_count = stage_clock()
        allocate(csf_subshell_occupations(max_orbitals, csf_count))
        allocate(csf_subshell_quantum_codes(max_orbitals,3,csf_count))
        allocate(csf_angular_codes(csf_angular_code_count, csf_count))
        call record_stage("csf_allocation", symmetry_block_index, stage_start_count)
        stage_start_count = stage_clock()
        call load_angular_block(input_options, core_subshell_count, symmetry_block_index)
        call record_stage("csf_load", symmetry_block_index, stage_start_count, int(csf_count,int64), storage=input_format_name)
        ! Open the temporary MCP channels for this symmetry block.
        call open_temporary_mcp_channels (symmetry_block_index, 30, files%temporary)
        ! Generate this symmetry block’s MCP data.
        call compute_mcp_angular_coefficients (symmetry_block_index, is_restart, mpi_rank, mpi_process_count, files%mcp, &
            & files%temporary)
    end do
    call close_angular_input(input_options)
    call record_stage("program_total_including_interaction_input", 0, program_start_count)
    call close_stage_profile()
    if (summary_output_unit /= -1) close(summary_output_unit)                                  ! Summary file
    close(angular_log_unit)                                  ! rangular.log
    if (debug .and. angular_debug_unit /= -1) close(angular_debug_unit)                       ! Debug file
    ! Report elapsed time.
    call finish_run_timer (start_clock_count, "RANGULAR")
    stop
end program rangular_main
