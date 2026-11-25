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
        use yes_no_input, only: read_yes_no
        use file_open, only: open_calculation_file

        ! Procedure interfaces
        use read_ci_restart_interface, only: read_ci_restart
        use read_ci_options_interface, only: read_ci_options
        implicit none (type, external)

        ! Arguments
        character  :: isotope_file_path*(*)
        character  :: radial_file_path*(*)
        character(len=8)  :: block_labels(*)

        logical, intent(out) :: is_restart

        ! Local constants
        character(len=11), parameter :: form = "UNFORMATTED"
        character(len=7), parameter :: status = "UNKNOWN"
        character(len=6), parameter :: restitle = "R92RES"

        ! Local declarations
        integer :: open_status, io_status
        logical :: found
        character :: r92res*6, file_path*11
        ! Compose the "rci.res" file name
        file_path = "rci.res"
        ! Ask if this is a restart
        if (nondefault_options_flag /= 0) then
            write (stderr_unit, *) "Restarting RCI90 ?"
            is_restart = read_yes_no()
        else
            is_restart = .false.
        end if
        ! Do some settings and checks
        if (is_restart) then
            ! ...Restart, make sure file exist
            inquire(file=file_path, exist=found)
            if (.not.found) stop "connect_ci_restart_file: .res does not exist"
        end if
        ! Open the .res file
        call open_calculation_file (restart_matrix_unit, file_path, form, status, open_status)
        if (open_status /= 0) stop "connect_ci_restart_file: Error openning .res file"
        ! If restart, load the contents. Otherwise generate them via getcid
        !
        ! But first of all, iccutblk() is needed in both cases
        call alloc (csf_cutoff_by_block, symmetry_block_count, "csf_cutoff_by_block", "connect_ci_restart_file")

        if (is_restart) then
            ! ...Check the signature of the file
            read (restart_matrix_unit, iostat=io_status) r92res
            if (io_status/=0 .or. r92res/=restitle) then
                close(restart_matrix_unit)
                stop "connect_ci_restart_file: Not RCI92 .res file"
            end if

            ! ...Read and check restart information
            call read_ci_restart

        else

            call read_ci_options (isotope_file_path, radial_file_path, block_labels)

        end if

        return
    end subroutine connect_ci_restart_file
end module connect_ci_restart_file_interface
