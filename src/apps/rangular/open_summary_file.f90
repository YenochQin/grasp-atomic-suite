!***********************************************************************
!                                                                      *
subroutine open_summary_file(file_path)
    ! Module dependencies
    use calculation_file_units, only: summary_output_unit
    !                                                                      *
    !   Open the  .sum  file on stream 24.                                 *
    !                                                                      *
    !   Call(s) to: [LIB92]:  OPENFL.                                      *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 11 Nov 1992   *
    !   Modified by Xinghong He               Last revision:  3 Jul 1998   *
    !
    !  File shared by mcpblk, mcpmpi
    !                                                                      *
    !***********************************************************************
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments
    character(len=*), intent(in) :: file_path

    ! Local declarations
    integer :: open_status
    character(len=120) :: filnam  ! The original did not compile
    character(len=11) :: form
    character(len=3)  :: status
    !-----------------------------------------------

    form = "FORMATTED"
    status = "NEW"
    !
    write (stderr_unit, *) "File ", file_path, " will be created as the", &
        " GENMCP SUMmary File;"
    write (stderr_unit, *) "enter another file name if this is not ", &
        "acceptable; null otherwise:"
    read (*, "(A)") filnam
    !
    if (len_trim(filnam) == 0) filnam = file_path
    !
    do
        call openfl (summary_output_unit, filnam, form, status, open_status)
        if (open_status == 0) exit
        do
            write (stderr_unit, *) "Enter a name for the GENMCP SUMmary", &
                " File that is to be created:"
            read (*, "(A)") filnam
            if (len_trim(filnam) /= 0) exit
        end do
    end do
    !
    return
end subroutine open_summary_file
