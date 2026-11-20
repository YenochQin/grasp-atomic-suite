module open_summary_file_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_summary_file
contains

    !***********************************************************************
    !                                                                      *
    subroutine open_summary_file(state_name, mixing_file_kind)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_file_units, only: hyperfine_output_unit, summary_output_unit
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
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: mixing_file_kind
        character, intent(in) :: state_name*24

        ! Local declarations
        integer   :: state_name_end, open_status
        character :: file_path*256, hyperfine_file_path*256, form*11, status*3
        !-----------------------------------------------
        !
        !   File  hfs92.sum  is FORMATTED
        !
        state_name_end = index(state_name," ")
        if (mixing_file_kind == 0) then
            file_path = state_name(1:state_name_end-1)//".ch"
            hyperfine_file_path = state_name(1:state_name_end-1)//".choffd"
        else
            file_path = state_name(1:state_name_end-1)//".h"
            hyperfine_file_path = state_name(1:state_name_end-1)//".hoffd"
        end if
        form = "FORMATTED"
        status = "NEW"
        !
        call open_calculation_file (hyperfine_output_unit, file_path, form, status, open_status)
        if (open_status /= 0) then
            write (output_unit, *) "Error when opening", file_path
            stop
        end if
        !
        call open_calculation_file (summary_output_unit, hyperfine_file_path, form, status, open_status)
        if (open_status /= 0) then
            write (output_unit, *) "Error when opening", hyperfine_file_path
            stop
        end if
        !
        return
    end subroutine open_summary_file
end module open_summary_file_interface
