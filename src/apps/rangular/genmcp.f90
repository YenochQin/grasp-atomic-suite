
!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**        *****   ******  **   **  **    **   *****   ******        ***
!**       **   **  **      ***  **  ***  ***  **   **  **   **       ***
!**       **       **      ***  **  ** ** **  **       **   **       ***
!**       **  ***  ****    ** ****  ** ** **  **       ******        ***
!**       **   **  **      **  ***  **    **  **       **            ***
!**       **   **  **      **   **  **    **  **   **  **            ***
!**        *****   ******  **   **  **    **   *****   **            ***
!**                                                                  ***
!**   Program for generating the energy expression for H(DC). This   ***
!**   program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,  ***
!**   and C. F. Fischer, 1990).                                      ***
!**                                                                  ***
!**                         GRASP92 Version                          ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      program genmcp
      use, intrinsic :: iso_fortran_env, only: output_unit
      use calculation_file_units, only: angular_log_unit, summary_output_unit
      use debug_output_settings, only: angular_debug_unit
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer                                    *
!                                                                      *
!   Modified by Gediminas Gaigalas:                         Feb 2017   *
!                               1) for new spin-angular integration,   *
!                               2) for sorting in the memory.          *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters,   only:  max_orbitals
      use iso_fortran_env, only: int64
      use rangular_input, only: angular_input_options, parse_angular_options, load_angular_header, &
                                load_angular_block, close_angular_input
      use rangular_profile, only: stage_clock, record_stage, close_stage_profile
      use array_memory_management, only: alloc
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use calculation_input_defaults, only: inpci, ndef, ndump
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calculation_state_files, only: state_files, state_name_length, files_for_state, select_state
      use run_timer_start, only: starttime
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use setsum_i, only: setsum

      use setmcp2_i, only: setmcp2
      use strsum_i, only: strsum
      use log_factorial_table, only: factt
!      USE settmp_I
      use settmpgg_i, only: settmpgg
      use run_timer_finish, only: stoptime

      use mcp_i, only: mcp
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nblk0 = 50   !-----------------------------------------------
      integer :: start_clock_count, core_subshell_count, symmetry_block_index
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical ::  debug, is_restart, yes
      character(len=128) :: startup_directory*128, permdir*128, tmpdir*128
      character(len=8):: idblk(nblk0)
      integer :: mpi_rank,mpi_process_count
      type(angular_input_options) :: input_options
      type(state_files) :: files
      character(state_name_length) :: state
      integer(int64) :: program_start_count, stage_start_count
      character(len=16) :: input_format_name
!-----------------------------------------------------------------------
!
      call parse_angular_options(input_options)
      call select_state(state, "Full interaction? (y/n)", yes)
      files = files_for_state(state)
      if (input_options%input_filename == "") input_options%input_filename = files%descriptors
      input_format_name = "descriptors"
      program_start_count = stage_clock()
      write(*,*)
      write(*,*) "RANGULAR"
      write(*,*) "This program performs angular integration "
      write(*,*) "Input file: ", input_options%input_filename
      write(*,*) "Output files: ", files%mcp, ".30, ", files%mcp, ".31, ..."
      write(*,*) files%angular_log
      write(*,*)
      open(newunit=angular_log_unit,file=files%angular_log,status="UNKNOWN")
      if (state /= "") write(angular_log_unit,"(A)") trim(state)
      mpi_rank = 0
      mpi_process_count = 1
      call starttime (start_clock_count, "RANGULAR")
!=======================================================================
!  Get NDEF
!=======================================================================
      if (mpi_rank == 0) then
!         WRITE (ISTDE, '(A)', ADVANCE='NO') 'Default settings?  (y/n) '
         if (yes) then
            ndef = 0
            write(angular_log_unit,"(A)") "y            ! Full interaction"
         else
            ndef = 1
            write(angular_log_unit,"(A)") "n            ! Full interaction"
         end if
      end if
!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    chkplt - compatibility of plant substitutions
!    setdbg - debug output control parameters
!    setmc - machine- and installation- dependent constants
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================
!CFF  delete chkplt
!     CALL CHKPLT ('GENMCP')
      call setdbg (debug, files%angular_debug)
      call setmc
!      IF (NDEF/=0 .AND. MYID==0) CALL SETSUM ('genmcp.sum')
      stage_start_count = stage_clock()
      call load_angular_header(input_options, core_subshell_count, nblk0, idblk)
      call record_stage("input_header", 0, stage_start_count, int(csf_count,int64), storage=input_format_name)
      is_restart = .false.
      stage_start_count = stage_clock()
      call setmcp2 (mpi_rank, mpi_process_count, core_subshell_count, idblk, files%mcp)
      if (ndef/=0 .and. mpi_rank==0) call strsum
      call factt
      call record_stage("setup_including_interaction_input", 0, stage_start_count)
!=======================================================================
!     For each block, generate and sort the data
!=======================================================================
      do symmetry_block_index = 1, symmetry_block_count
         csf_count = csf_counts_by_block(symmetry_block_index)                        ! This ncf goes to common
         if (mpi_rank == 0) then
            write (output_unit, *)
            write (output_unit, *) "Block ", symmetry_block_index, ",  ncf = ", csf_count
         end if
            !*** Load current CSL block. Memories de-allocated in mcp ***
         stage_start_count = stage_clock()
         call alloc (csf_subshell_occupations, max_orbitals, csf_count, "IQA", "GENMCP")
         call alloc (csf_subshell_quantum_codes, max_orbitals,3,csf_count, "JQSA", "GENMCP")
         call alloc (csf_coupling_quantum_codes, max_orbitals, csf_count, "JCUPA", "GENMCP")
         call record_stage("csf_allocation", symmetry_block_index, stage_start_count)
!
         stage_start_count = stage_clock()
         call load_angular_block(input_options, core_subshell_count, symmetry_block_index)
         call record_stage("csf_load", symmetry_block_index, stage_start_count, int(csf_count,int64), storage=input_format_name)
         !*** Open tmp.xx files for block nb ***
!GG         CALL SETTMP (NB, KMAX, 'tmp')
         call settmpgg (symmetry_block_index, 30, files%temporary)
         !*** Generation of MCP coefficients ***
         call mcp (symmetry_block_index, is_restart, mpi_rank, mpi_process_count, files%mcp, files%temporary)
      end do
      call close_angular_input(input_options)
      call record_stage("program_total_including_interaction_input", 0, program_start_count)
      call close_stage_profile()
      if (summary_output_unit /= -1) close(summary_output_unit)                                  ! Summary file
      close(angular_log_unit)                                  ! rangular.log
      if (debug .and. angular_debug_unit /= -1) close(angular_debug_unit)                       ! Debug file
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      call stoptime (start_clock_count, "RANGULAR")
      stop
      end program genmcp
