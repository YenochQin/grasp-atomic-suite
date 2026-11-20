module open_summary_file_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_summary_file
contains

    !***********************************************************************
    !                                                                      *
    subroutine open_summary_file(name, nci)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_file_units, only: summary_output_unit
        !                                                                      *
        !   Open the  .sum  files on stream 24 and 29                          *
        !                                                                      *
        !   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
        !                                                                      *
        !   Updated by Per Jonsson                               28 Oct 1999   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Gediminas Gaigalas 11/18/19
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nci
        character, intent(in) :: name*24

        ! Local declarations
        integer   :: k, open_status
        character :: file_path*256, form*11, status*3
        !-----------------------------------------------
        !
        !   File  rdensity.sum  is FORMATTED
        !
        k = index(name," ")
        if (nci == 0) then
            file_path = name(1:k-1)//".ci"
        else
            file_path = name(1:k-1)//".i"
        end if
        form = "FORMATTED"
        status = "NEW"
        !
        call open_calculation_file (summary_output_unit, file_path, form, status, open_status)
        if (open_status /= 0) then
            write (output_unit, *) "Error when opening", file_path
            stop
        end if
        !
        return
    end subroutine open_summary_file
end module open_summary_file_interface
