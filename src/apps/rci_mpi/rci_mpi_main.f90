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
program rci_mpi_main
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
    use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_process_count, mpi_rank, processor_name, &
        processor_name_length
    use hamiltonian_correction_flags, only: include_transverse_interaction
    use calculation_state_files, only: read_state_name, orbital_input_path, descriptor_input_path
    use yes_no_input, only: read_yes_no
    use mpi_csf_block_input, only: open_csf_blocks_mpi

    ! Procedure interfaces
    use configure_ci_debug_output_interface, only: configure_ci_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file
    use connect_ci_restart_file_interface, only: connect_ci_restart_file
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
    use mpi_utilities, only: finalize_mpi_runtime_with_status
    implicit none (type, external)
    integer, parameter :: maximum_input_symmetry_blocks = 50
    integer, parameter :: nblk=50
    integer :: start_clock_count
    integer :: core_subshell_count, max_orbital_two_j
    logical :: yes
    !      CHARACTER(LEN=128) :: NAME, TMPDIR, PERMDIR, ISOFILE
    character(len=128) :: startup_directory, state_name, scratch_directory, permanent_directory
    character(:), allocatable :: name, isotope_file_path, orbital_path, descriptor_file_path
    character(len=8) :: block_labels(maximum_input_symmetry_blocks)
    !-----------------------------------------------
    !=======================================================================
    !  Start mpi --- get processor info: myid, nprocs, host name; and print
    !=======================================================================

    call initialize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length, &
        & start_clock_count,       &
        startup_directory, permanent_directory, scratch_directory, "RCI_MPI")


    !=======================================================================
    !  Get NDEF on node-0 and then send to all nodes
    !=======================================================================

    if (mpi_rank == 0) then
        write (stderr_unit,*)
        write (stderr_unit, *) "Default settings? "
        yes = read_yes_no()
        if (yes) then
            nondefault_options_flag = 0
        !cjb fort.734
        !cjb        write(734,'(A)') 'y            ! Default settings'
        else
            nondefault_options_flag = 1
            !cjb fort.734
            !cjb        write(734,'(A)') 'n            ! Default settings'
        end if
    end if
    call mpi_bcast (nondefault_options_flag,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    !     CALL MPI_Bcast (imethod,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


    !=======================================================================
    !  Get name of the state (used in files like <name>.c, <name>.s)
    !  info used only on node-0
    !=======================================================================

    name = ""
    isotope_file_path = ""
    orbital_path = ""
    descriptor_file_path = ""
    if (mpi_rank == 0) then
        call read_state_name(state_name)
        open(newunit=calculation_log_unit,file=trim(state_name)//".clog",status="unknown")
        write(calculation_log_unit,"(a)") "y            ! Default settings"
        write(calculation_log_unit,"(a)") trim(state_name)
        isotope_file_path = trim(permanent_directory)//"/isodata"
        name = trim(permanent_directory)//"/"//trim(state_name)
        orbital_path = orbital_input_path(name)
        descriptor_file_path = descriptor_input_path(name)
    end if

    !=======================================================================
    do


        prerun_stage = 0
        !=======================================================================

        !=======================================================================
        !
        !  Checks and settings... Mostly done in backyard.
        !
        !    SETDBG - Debug output control parameters - all set to false
        !    SETMC - machine- and installation- dependent constants
        !    SETCON - physical constants
        !    SETSUM - open the summary file
        !    cslhmpi - load header of the CSL file
        !    SETRES - setup the .res  file, one for each node
        !    SETMIXmpi - set up mixing coefficients file AND process user input
        !    FACTT - table of logarithms of factorials setup
        !    SETCSLmpi - open, check, load data from, and close the  .csl  file
        !=======================================================================

        call configure_ci_debug_output
        call initialize_machine_constants
        call initialize_physical_constants
        if (mpi_rank == 0)  call open_summary_file (name)
        call open_csf_blocks_mpi (descriptor_file_path, core_subshell_count, block_labels)
        call connect_ci_restart_file (isotope_file_path, orbital_path, block_labels)
        call open_mixing_output_mpi (name, block_labels)
        if (mpi_rank == 0)  call write_input_summary
        call initialize_log_factorial_table

        !=======================================================================
        !   Calculate all the needed Rk integrals - genintrkwrap calls genintrk
        !   to compute the integrals (each processor allocates the same amount
        !   of memory, but only computes part of the integrals); and then do the
        !   "MPI_Allreduce" stuff.
        !=======================================================================

        call tabulate_slater_integrals_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)

        !
        !   If transverse interaction comput Breit integrals of type 1 and 2
        !
        if (include_transverse_interaction) then
            call tabulate_breit_integrals_1_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)
            call tabulate_breit_integrals_2_mpi (mpi_rank, mpi_process_count, max_orbital_two_j)
        end if

        !=======================================================================
        !   Proceed with the CI calculation
        !=======================================================================
        write (output_unit, *) "Calling MATRIX..."

        call solve_hamiltonian_blocks (core_subshell_count, max_orbital_two_j)

        if (prerun_stage /= 1) exit
        prerun_stage = 2
    end do
    !=======================================================================
    !  Execution finished; Statistics output
    !=======================================================================

    call finalize_mpi_runtime_with_status (mpi_rank, mpi_process_count, processor_name, processor_name_length,          &
        start_clock_count, "RCI_MPI")

    stop
end program rci_mpi_main
