! Open the mixing output file and write its program-specific header.
module open_mixing_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_mixing_output
contains

    ! Opens the  .mix  file on stream 25; writes a header to this file.
    !
    ! Written by Farid A. Parpia. Last revision: 24 Dec 1992
    ! Modified by Xinghong He. Last revision: 13 Jul 1998
    subroutine open_mixing_output(name)
        ! Module dependencies
        use calculation_file_units, only: mixing_output_unit
        use mixing_file_format, only: open_mixing_stream
        use mpi_io_errors, only: fail_input_io_mpi
        implicit none (type, external)

        ! Arguments
        character, intent(in) :: name*(*)

        ! Local declarations
        integer :: open_status
        call open_mixing_stream(name, mixing_output_unit, open_status)
        if (open_status /= 0) call fail_input_io_mpi(open_status, mixing_output_unit, &
            "Cannot create RMCDHF mixing stream: "//trim(name))
        return
    end subroutine open_mixing_output
end module open_mixing_output_interface
