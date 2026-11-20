module configure_debug_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_debug_output
contains

    !***********************************************************************
    !                                                                      *
    subroutine configure_debug_output(default_file_path)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit
        !                                                                      *
        !   This subroutine sets the arrays that control debug printout from   *
        !   the radial and angular modules of the GRASP92 suite.               *
        !                                                                      *
        !   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
        !                                                                      *
        !   Written by Farid A Parpia               Last update: 15 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use debug_output_settings, only: ldbpa, ldbpg, rmcdhf_debug_flags
        use calculation_input_defaults, only: nondefault_options_flag
        use standard_io_units, only: stderr_unit
        use yes_no_input, only: read_yes_no
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in) :: default_file_path

        ! Local declarations
        integer :: open_status
        logical :: yes
        character :: file_path*256, form*11, status*3
        !-----------------------------------------------
        !
        !
        !   Initialise the arrays that control the debug printout
        !
        ldbpa = .false.
        !
        ldbpg = .false.
        !
        rmcdhf_debug_flags = .false.
        !
        if (nondefault_options_flag == 0) return

        write (stderr_unit, *) "Generate debug printout?"
        yes = read_yes_no()
        if (yes) then
            !
            !   The  .dbg  file is formatted; open it on unit 99
            !
            form = "FORMATTED"
            status = "NEW"
            !
            write (stderr_unit, *) "File ", default_file_path, " will be created as the ", &
                "ERWF DeBuG Printout File; "
            write (stderr_unit, *) "enter another file name if this is not ", &
                "acceptable; null otherwise:"
            read (*, "(A)") file_path
            !
            if (len_trim(file_path) == 0) file_path = default_file_path
            !
            do
                call open_calculation_file (angular_debug_unit, file_path, form, status, open_status)
                if (open_status == 0) exit
                do
                    write (stderr_unit, *) "Enter a name for the ERWF DeBuG Printout ", &
                        "file that is to be created:"
                    read (*, "(A)") file_path
                    if (len_trim(file_path) > 0) exit
                end do
            end do
            !
            !   Set options for general printout
            !
            write (stderr_unit, *) " Print out the machine constants used?"
            yes = read_yes_no()
            if (yes) ldbpg(1) = .true.
            write (stderr_unit, *) " Print out the physical constants used?"
            yes = read_yes_no()
            if (yes) ldbpg(2) = .true.
            !
            !   Set options for radial modules
            !
            write (stderr_unit, *) " Printout from RADGRD?"
            yes = read_yes_no()
            if (yes) rmcdhf_debug_flags(1) = .true.
            write (stderr_unit, *) " Printout from NUCPOT?"
            yes = read_yes_no()
            if (yes) rmcdhf_debug_flags(2) = .true.
            write (stderr_unit, *) " Printout from TFPOT?"
            yes = read_yes_no()
            if (yes) rmcdhf_debug_flags(26) = .true.
            !
            !   Set options for angular modules
            !
            write (stderr_unit, *) " Printout from LODCSL?"
            yes = read_yes_no()
            if (yes) ldbpa(1) = .true.
            !
        end if
        !
        return
    end subroutine configure_debug_output
end module configure_debug_output_interface
