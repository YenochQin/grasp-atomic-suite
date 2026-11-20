

!***********************************************************************
!                                                                      *
subroutine open_summary_file(default_file_path)
    ! Module dependencies
    use calculation_file_units, only: summary_output_unit
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    !                                                                      *
    !   Open the  .sum  file on stream 24.                                 *
    !                                                                      *
    !   Call(s) to: [LIB92]:  OPENFL.                                      *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    !
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments
    character(len=*), intent(in) :: default_file_path

    ! Local declarations
    integer :: open_status
    character :: file_path*256, form*11, status*3
    !-----------------------------------------------
    !
    !   File  erwf.sum  is FORMATTED
    !
    form = "FORMATTED"
    status = "NEW"
    !
    write (stderr_unit, *) "File ", default_file_path, " will be created as the ", &
        "ERWF SUMmary File; "
    write (stderr_unit, *) "enter another file name if this is not ", &
        "acceptable; null otherwise:"
    read (*, "(A)") file_path
    !
    if (len_trim(file_path) == 0) file_path = default_file_path
    !
    do
        call openfl (summary_output_unit, file_path, form, status, open_status)
        if (open_status == 0) exit
        do
            write (stderr_unit, *) "Enter a name for the ERWF SUMmary File ", &
                "that is to be created:"
            read (*, "(A)") file_path
            if (len_trim(file_path) > 0) exit
        end do
    end do
    !
    return
end subroutine open_summary_file
