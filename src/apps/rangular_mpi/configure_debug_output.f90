module configure_debug_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_debug_output
contains

    !***********************************************************************
    !                                                                      *
    subroutine configure_debug_output (debug, file_path)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit
        !                                                                      *
        !   This routine opens the  .dbg  file and sets the arrays that con-   *
        !   trol debug printout from the GENMCP program.                       *
        !                                                                      *
        !   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
        !                                                                      *
        !   Written by Farid A Parpia               Last update: 08 Dec 1992   *
        !   Modified by Xinghong He                 Last update: 03 Jul 1998   *
        !
        !   Used by mcpvu, mcpmpivu
        !                                                                      *
        !***********************************************************************
        use debug_output_settings, only: ldbpa, ldbpg
        use calculation_input_defaults, only: nondefault_options_flag
        use standard_io_units, only: stderr_unit
        use yes_no_input, only: read_yes_no
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        logical :: debug
        character(len=*) :: file_path

        ! Local declarations
        integer :: i, open_status
        !
        character(len = len (file_path)) :: filnam
        character(len=11) :: form
        character(len=3) :: status
        !
        !   Initialise the arrays that control the debug printout
        !
        do i = 1, 5
            ldbpa(i) = .false.
            ldbpg(i) = .false.
        end do
        !
        return      ! Skip these questions and return immediately
        !
        if (nondefault_options_flag == 0) then
            return
        end if

        write (stderr_unit,*) "Generate debug printout?"
        debug = read_yes_no ()
        if (debug) then
            !
            !   The  .dbg  file is formatted; open it on unit 99
            !
            form = "FORMATTED"
            status = "NEW"
            !
            write (stderr_unit,*) "File ", file_path," will be created as the"    &
                &,                 " GENMCP DeBuG Printout File; "
            write (stderr_unit,*) " enter another file name if this is not "     &
                &,                 "acceptable; null otherwise:"
            read (*,"(A)") filnam
            !
            if (len_trim (filnam) == 0) filnam = file_path
            !
            do
                call open_calculation_file (angular_debug_unit, filnam, form, status, open_status)
                if (open_status == 0) exit
                do
                    write (stderr_unit,*) "Enter a name for the GENMCP DeBuG Printout"&
                        &,                    " file that is to be created:"
                    read (*,"(A)") filnam
                    if (len_trim(filnam) /= 0) exit
                end do
            end do
            !
            !   Set options for general printout
            !
            write (stderr_unit,*) " Print out the machine constants used?"
            ldbpg(1) = read_yes_no ()
            !
            !   Set options for angular modules
            !
            write (stderr_unit,*) " Printout from LODCSL? (Not used)"
            ldbpa(1) = read_yes_no ()
            write (stderr_unit,*) " Print out T coefficients?"
            ldbpa(2) = read_yes_no ()
            write (stderr_unit,*) " Print out Coulomb V coefficients?"
            ldbpa(3) = read_yes_no ()
            write (stderr_unit,*) " Print out sparse matrix definition arrays?"
            ldbpa(4) = read_yes_no ()

        end if
        !
        return
    end subroutine configure_debug_output
end module configure_debug_output_interface
