!***********************************************************************
!                                                                      *
subroutine open_summary_file(name, nci)
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
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: nci
    character, intent(in) :: name*24

    ! Local declarations
    integer   :: k, open_status
    character :: file_path*256, filnam2*256, default_file_path*11, form*11, status*3
    !-----------------------------------------------
    !
    !   File  hfs92.sum  is FORMATTED
    !
    k = index(name," ")
    if (nci == 0) then
        file_path = name(1:k-1)//".ch"
        filnam2 = name(1:k-1)//".choffd"
    else
        file_path = name(1:k-1)//".h"
        filnam2 = name(1:k-1)//".hoffd"
    end if
    form = "FORMATTED"
    status = "NEW"
    !
    call openfl (hyperfine_output_unit, file_path, form, status, open_status)
    if (open_status /= 0) then
        write (output_unit, *) "Error when opening", file_path
        stop
    end if
    !
    call openfl (summary_output_unit, filnam2, form, status, open_status)
    if (open_status /= 0) then
        write (output_unit, *) "Error when opening", filnam2
        stop
    end if
    !
    return
end subroutine open_summary_file
