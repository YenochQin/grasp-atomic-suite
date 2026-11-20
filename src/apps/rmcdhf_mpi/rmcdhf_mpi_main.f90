!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**       ******    *****    *****   *******   *****    *****        ***
!**       **   **  **   **  **   **  **       **   **  **   **       ***
!**       **   **  **       **       **       **   **       **       ***
!**       ******    *****   **       ****      *****       **        ***
!**       **  **        **  **       **          **       **         ***
!**       **   **  **   **  **   **  **         **      **           ***
!**       **   **   *****    *****   **        **      *******       ***
!**                                                                  ***
!**            Relativistic Self-Consistent-Field Program            ***
!**                                                                  ***
!**   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
!**   Grant, and C. F. Fischer, 1990)                                ***
!**                                                                  ***
!**                            GRASP92                               ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!**   Modified by C.F. FIscher for block input                       ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
program rmcdhf_mpi_main
    ! Module dependencies
    use calculation_file_units, only: calculation_log_unit
    !                                                                      *
    !   Entry routine for RSCFVU. Controls the entire computation.      *
    !                                                                      *
    !   Call(s) to: [LIB92]: SETMC, SETCON.                                *
    !               [RSCF92]: CHKPLT, setcsl, setdbg
    !                        SETMIX, SETRES, SETSUM, STRSUM.
    !               [NJGRAF]: FACTT.                                       *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
    !   Block version by Xinghong He          Last revision: 17 Aug 1998   *
    !   Midified by G. Gaigalas                              05 Feb 2017   *
    !      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
    !                                  JCUPA(NNNWP*NCF)                    *
    !                                                                      *
    !***********************************************************************
    use calculation_input_defaults, only: inpci, ndef, ndump
    use closed_core_subshells, only: closed_core_subshell_count
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use mpi_state, only: mpi_process_count, mpi_rank, &
        & processor_name, processor_name_length
    !GG po to isimti
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use calculation_state_files, only: check_orbital_output
    use calculation_state_files, only: state_files, state_name_length, files_for_state
    use mpi_calculation_state, only: select_mpi_state

    ! Procedure interfaces
    use configure_debug_output_mpi_interface, only: configure_debug_output_mpi

    ! Module dependencies
    use machine_constants, only: setmc
    use physical_constants, only: setcon

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file
    use connect_mcp_files_interface, only: connect_mcp_files
    use load_csf_basis_data_mpi_interface, only: load_csf_basis_data_mpi
    use read_scf_options_mpi_interface, only: read_scf_options_mpi
    use write_input_summary_interface, only: write_input_summary
    use open_mixing_output_interface, only: open_mixing_output

    ! Module dependencies
    use log_factorial_table, only: factt

    ! Procedure interfaces
    use run_scf_iterations_mpi_interface, only: run_scf_iterations_mpi

    ! Module dependencies
    use mpi_utilities, only: startmpi2
    use mpi_utilities, only: stopmpi2
    implicit none (type, external)
    integer, parameter :: nblk0 = 50  ! Maximum number of blocks
    integer :: ncore1, start_clock_count, lenperm
    type(state_files) :: files
    character(state_name_length) :: state
    logical :: extended_optimal_level, yes
    character :: block_labels(nblk0)*8
    !      CHARACTER (LEN = *) :: host
    character (len = 128) :: startup_directory,serial_io_directory,rank_scratch_root
    character(:), allocatable :: file_rcsl, file1, file2
    !-----------------------------------------------
    !

    !=======================================================================
    !  Start mpi --- get processor info: myid, nprocs, host name; and print
    !=======================================================================
    startup_directory = "  "  ;    file_rcsl = "  "
    serial_io_directory = "  "   ;    file1     = "  "
    rank_scratch_root = "  "    ;    file2     = "  "
    call startmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count,       &
        startup_directory, serial_io_directory, rank_scratch_root, "RMCDHF_MPI")
    lenperm = len_trim (serial_io_directory)

    !=======================================================================
    !  Get NDEF
    !=======================================================================

    call select_mpi_state(state, "Default settings? (y/n)", yes, serial_io_directory)
    files = files_for_state(state)
    if (mpi_rank == 0) call check_orbital_output(trim(serial_io_directory)//"/"//files%orbitals, writing_initial=.false.)
    ndef = merge(0, 1, yes)
    if (mpi_rank == 0) then
        open(newunit=calculation_log_unit,file=trim(startup_directory)//"/"//files%scf_log,status="UNKNOWN")
        if (state /= "") then
            write(calculation_log_unit,"(A)") trim(state)
            write(calculation_log_unit,"(A)") merge("y","n",yes)//"            ! Default settings"
        end if
    end if

    !=======================================================================
    !
    !  Checks and settings... Mostly done in backyard.
    !
    !    CHKPLT - compatibility of plant substitutions
    !    SETDBGmpi - Debug output control parameters
    !    SETMC - machine- and installation- dependent constants
    !    SETCON - physical constants
    !    SETSUM - open the summary file
    !    SETMCP - open and check the  .mcp  files
    !    SETCSLmpi - open, check, load data from, and close the  .csl  file
    !    STRSUM - append a summary of the inputs to the  .sum  file
    !    SETMIX - mixing coefficients file setup
    !    FACTT - table of logarithms of factorials setup
    !=======================================================================
    call configure_debug_output_mpi (trim(serial_io_directory)//"/"//files%scf_debug)
    call setmc
    call setcon

    if (mpi_rank == 0) call open_summary_file (trim(serial_io_directory)//"/"//files%scf_summary)

    call connect_mcp_files (closed_core_subshell_count, nblk0, block_labels, files%mcp)
    if(mpi_rank == 0) file_rcsl = trim(serial_io_directory)//"/"//files%descriptors
    call load_csf_basis_data_mpi (file_rcsl, ncore1, block_labels)
    if (closed_core_subshell_count /= ncore1) error stop "rscfmpivu: ncore"

    !=======================================================================
    !  Gather all remaining information and perform some setup. This
    !  part (routine) asks for user-inputs.
    !=======================================================================

    if(mpi_rank == 0) file1 = serial_io_directory(1:lenperm) // "/isodata"
    if(mpi_rank == 0) file2 = trim(serial_io_directory)//"/"//files%initial_orbitals
    call read_scf_options_mpi (extended_optimal_level, block_labels, file1, file2 )
    file1 = "  "
    if(mpi_rank == 0) file1 = trim(serial_io_directory)//"/"//files%mixing
    if (mpi_rank == 0) then
        call write_input_summary
        if (extended_optimal_level) call open_mixing_output (file1)
    end if

    call factt

    !=======================================================================
    !  Proceed with the SCF calculation close all files except
    !  the  .sum  file
    !=======================================================================

    file1 = "  "
    if(mpi_rank == 0) file1 = trim(serial_io_directory)//"/"//files%orbitals

    call run_scf_iterations_mpi (extended_optimal_level, file1)

    !=======================================================================
    !  Execution finished; Statistics output
    !=======================================================================

    call stopmpi2 (mpi_rank, mpi_process_count, processor_name, processor_name_length, start_clock_count, "RMCDHF_MPI")

end program rmcdhf_mpi_main
