! Open the shared mixing output on the designated MPI rank.
module open_mixing_output_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_mixing_output_mpi
contains

    ! Opens the  .mix  file on stream 25; writes a header to this file;
    ! calls LODMIX to interactively determine the eigenpairs required.
    !
    ! Written by Farid A. Parpia. Last revision: 18 Dec 1992
    ! Modified by Xinghong He. Last revision: 23 Jun 1998
    subroutine open_mixing_output_mpi(name, block_labels)
        ! Module dependencies
        use calculation_state_files, only: ci_mixing_suffix
        use calculation_file_units, only: mixing_output_unit
        use mpi_state, only: mpi_error_code, mpi_rank
        use file_open, only: open_calculation_file

        ! Procedure interfaces
        use read_ci_mixing_coefficients_mpi_interface, only: read_ci_mixing_coefficients_mpi

        ! Module dependencies
        use mpi_utilities, only: finalize_mpi_runtime
        implicit none (type, external)

        ! Arguments
        character(len=*) , intent(in) :: name
        character(len=8) :: block_labels(*)

        ! Local constants
        character(len=11), parameter :: form = "UNFORMATTED"
        character(len=7), parameter  :: status = "UNKNOWN"

        if (mpi_rank == 0) then
            call open_calculation_file (mixing_output_unit, trim(name)//ci_mixing_suffix, form, status, mpi_error_code)
            if (mpi_error_code /= 0) then
                call finalize_mpi_runtime ("setmix: Error when opening .cimix file", mpi_rank)
            end if

            write (mixing_output_unit) "G92MIX"
        end if

        call read_ci_mixing_coefficients_mpi (block_labels)

        return
    end subroutine open_mixing_output_mpi
end module open_mixing_output_mpi_interface
