module open_summary_file_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_summary_file
contains

    !    Open the  .sum  files on stream 24 and 29
    !
    ! Written by Farid A. Parpia. Last revision: 24 Dec 1992
    !
    !    Updated by Per Jonsson                               28 Oct 1999
    !
    ! Translated by Gediminas Gaigalas 11/18/19
    subroutine open_summary_file(name, nci)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_state_files, only: isotope_shift_report_path
        use calculation_file_units, only: summary_output_unit
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nci
        character, intent(in) :: name*24

        ! Local declarations
        integer :: open_status
        character(:), allocatable :: file_path

        file_path = isotope_shift_report_path(name, nci == 0)
        call open_calculation_file (summary_output_unit, file_path, "FORMATTED", "NEW", open_status)
        if (open_status /= 0) then
            write (output_unit, *) "Error when opening", file_path
            stop
        end if
        return
    end subroutine open_summary_file
end module open_summary_file_interface
