! Open and validate the CI restart file, or initialize its records for a new run.
module connect_ci_restart_file_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: connect_ci_restart_file
contains

    !   Open, check, load data from the  .res  file.
    !
    ! Written by Farid A. Parpia. Last revision: 06 Oct 1992
    ! Modified by Xinghong. Last revision: 23 Jun 1998
    subroutine connect_ci_restart_file(isotope_file_path, radial_file_path, block_labels, is_restart)
        ! Module dependencies
        use array_memory_management, only: alloc
        use calculation_input_defaults, only: nondefault_options_flag
        use ci_restart_units, only: restart_matrix_unit
        use hamiltonian_symmetry_blocks, only: symmetry_block_count
        use csf_block_cutoffs, only: csf_cutoff_by_block
        use standard_io_units, only: stderr_unit
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank
        use yes_no_input, only: read_yes_no
        use rci_restart_records, only: open_rci_restart_file

        ! Procedure interfaces
        use read_ci_restart_interface, only: read_ci_restart
        use read_ci_options_interface, only: read_ci_options

        use mpi_io_errors, only: abort_input_mpi
        implicit none (type, external)

        ! Arguments
        character  :: isotope_file_path*(*)
        character  :: radial_file_path*(*)
        character(len=8)  :: block_labels(*)

        logical, intent(out) :: is_restart

        character(11) :: file_path
        character(3) :: block_label_text
        ! Compose the "rciXXX.res" file names (for each node)
        write (block_label_text, "(I3.3)") mpi_rank
        file_path = "rci" // block_label_text // ".res"
        ! Ask if this is a restart
        if (mpi_rank == 0) then
            if (nondefault_options_flag /= 0) then
                write (stderr_unit, *) "Restarting RCI90 ?"
                is_restart = read_yes_no()
            else
                is_restart = .false.
            end if

        end if
        call mpi_bcast (is_restart, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        call open_rci_restart_file(restart_matrix_unit, file_path, is_restart, abort_input_mpi)
        ! If restart, load the contents. Otherwise generate them via getcid
        !
        ! But first of all, iccutblk() is needed in both cases
        call alloc (csf_cutoff_by_block, symmetry_block_count, "csf_cutoff_by_block", "connect_ci_restart_file")

        if (is_restart) then
            call read_ci_restart

        else

            call read_ci_options (isotope_file_path, radial_file_path, block_labels)

        end if

        return
    end subroutine connect_ci_restart_file
end module connect_ci_restart_file_interface
