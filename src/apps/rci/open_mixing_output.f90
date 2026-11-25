! Open the mixing output file and write its program-specific header.
module open_mixing_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_mixing_output
contains

    ! Opens the  .mix  file on stream 25; writes a header to this file;
    ! calls LODMIX to interactively determine the eigenpairs required.
    !
    ! Written by Farid A. Parpia. Last revision: 18 Dec 1992
    ! Modified by Xinghong He. Last revision: 23 Jun 1998
    subroutine open_mixing_output(name, block_labels)
        ! Module dependencies
        use calculation_state_files, only: ci_mixing_suffix
        use calculation_file_units, only: mixing_output_unit
        use file_open, only: open_calculation_file

        ! Procedure interfaces
        use read_ci_mixing_coefficients_interface, only: read_ci_mixing_coefficients
        implicit none (type, external)

        ! Arguments
        character(len=*) , intent(in) :: name
        character(len=8), intent(in) :: block_labels(*)
        integer :: mixing_header(6)

        ! Local constants
        character(len=11), parameter :: form = "UNFORMATTED"
        character(len=7), parameter :: status = "UNKNOWN"

        ! Local declarations
        integer :: k, open_status

        call read_ci_mixing_coefficients (block_labels, mixing_header)

        k = index(name," ")
        call open_calculation_file (mixing_output_unit, name(1:k-1)//ci_mixing_suffix, form, status, open_status)
        if (open_status /= 0) stop "open_mixing_output: Error when opening .cimix file"

        write (mixing_output_unit) "G92MIX"
        write (mixing_output_unit) mixing_header

        return
    end subroutine open_mixing_output
end module open_mixing_output_interface
