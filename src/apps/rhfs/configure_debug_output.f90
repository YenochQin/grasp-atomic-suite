module configure_debug_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_debug_output
contains

    !***********************************************************************
    !                                                                      *
    subroutine configure_debug_output
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit
        !-----------------------------------------------
        !                                                                      *
        !   This subroutine sets the arrays that control debug printout from   *
        !   the radial and angular modules of the GRASP92 suite.               *
        !                                                                      *
        !   Call(s) to: [LIB92]: GETYN, LENGTH, OPENFL.                        *
        !                                                                      *
        !   Written by Farid A Parpia               Last update: 24 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use debug_output_settings, only: ldbpa, ldbpg, rmcdhf_debug_flags
        use calculation_input_defaults, only: nondefault_options_flag
        use yes_no_input, only: read_yes_no
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Local declarations
        integer   :: open_status
        logical   :: prompt_answer
        character :: file_path*256, default_file_path*11, form*11, status*3
        !-----------------------------------------------
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

        write (output_unit, *) "Generate debug printout?"
        prompt_answer = read_yes_no()
        if (prompt_answer) then
            !
            !   The  .dbg  file is formatted; open it on unit 99
            !
            default_file_path = "hfs92.dbg"
            form = "FORMATTED"
            status = "NEW"
            !
            write (output_unit, *) "File  hfs92.dbg  will be created as the"
            write (output_unit, *) " HFS92 DeBuG Printout File; enter another"
            write (output_unit, *) " file name if this is not acceptable;"
            write (output_unit, *) " null otherwise:"
            read (*, "(A)") file_path
            !
            if (len_trim(file_path) == 0) file_path = default_file_path
            !
            do
                call open_calculation_file (angular_debug_unit, file_path, form, status, open_status)
                if (open_status == 0) exit
                do
                    write (output_unit, *) "Enter a name for the HFS92 DeBuG Printout"
                    write (output_unit, *) " file that is to be created:"
                    read (*, "(A)") file_path
                    if (len_trim(file_path) /= 0) exit
                end do
            end do
            !
            !   Set options for general printout
            !
            write (output_unit, *) " Print out the machine constants used?"
            prompt_answer = read_yes_no()
            if (prompt_answer) ldbpg(1) = .true.
            write (output_unit, *) " Print out the physical constants used?"
            prompt_answer = read_yes_no()
            if (prompt_answer) ldbpg(2) = .true.
            !
            !   Set options for radial modules
            !
            write (output_unit, *) " Printout from radial modules?"
            prompt_answer = read_yes_no()
            if (prompt_answer) then
                write (output_unit, *) " Printout from RADGRD?"
                prompt_answer = read_yes_no()
                if (prompt_answer) rmcdhf_debug_flags(1) = .true.
                write (output_unit, *) " Printout from NUCPOT?"
                prompt_answer = read_yes_no()
                if (prompt_answer) rmcdhf_debug_flags(2) = .true.
                write (output_unit, *) " Printout from LODRWF?"
                prompt_answer = read_yes_no()
                if (prompt_answer) rmcdhf_debug_flags(3) = .true.
                !
            end if
            !
            !   Set options for angular modules
            !
            write (output_unit, *) " Printout from angular modules?"
            prompt_answer = read_yes_no()
            if (prompt_answer) then
                write (output_unit, *) " Printout from LODCSL?"
                prompt_answer = read_yes_no()
                if (prompt_answer) ldbpa(1) = .true.
                write (output_unit, *) " Print out T coefficients?"
                prompt_answer = read_yes_no()
                if (prompt_answer) ldbpa(2) = .true.
            end if
            !
        end if
        !
        return
    end subroutine configure_debug_output
end module configure_debug_output_interface
