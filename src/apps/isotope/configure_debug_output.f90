! Set diagnostic output options for the calculation.
module configure_debug_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_debug_output
contains

    !    This subroutine sets the arrays that control debug printout from
    !    the radial and angular modules of the GRASP92 suite.
    !
    ! Written by Farid A Parpia. Last update: 24 Dec 1992
    !
    ! Translated by Gediminas Gaigalas 11/18/19
    subroutine configure_debug_output
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit
        use debug_output_settings, only: ldbpa, ldbpg, rmcdhf_debug_flags
        use calculation_input_defaults, only: nondefault_options_flag
        use yes_no_input, only: read_yes_no
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Local declarations
        integer   :: open_status
        logical   :: yes
        character :: file_path*256, default_file_path*11, form*11, status*3
        ! Initialize the diagnostic output flags.
        ldbpa = .false.
        ldbpg = .false.
        rmcdhf_debug_flags = .false.
        if (nondefault_options_flag == 0) return

        write (output_unit, *) "Generate debug printout?"
        yes = read_yes_no()
        if (yes) then
            ! The  .dbg  file is formatted; open it on unit 99
            default_file_path = "isotope.dbg"
            form = "FORMATTED"
            status = "NEW"
            write (output_unit, *) "File  isotope.dbg  will be created as the"
            write (output_unit, *) " isotope DeBuG Printout File; enter another"
            write (output_unit, *) " file name if this is not acceptable;"
            write (output_unit, *) " null otherwise:"
            read (*, "(A)") file_path
            if (len_trim(file_path) == 0) file_path = default_file_path
            do
                call open_calculation_file (angular_debug_unit, file_path, form, status, open_status)
                if (open_status == 0) exit
                do
                    write (output_unit, *) "Enter a name for the isotope DeBuG Printout"
                    write (output_unit, *) " file that is to be created:"
                    read (*, "(A)") file_path
                    if (len_trim(file_path) /= 0) exit
                end do
            end do
            ! Read the general output options.
            write (output_unit, *) " Print out the machine constants used?"
            yes = read_yes_no()
            if (yes) ldbpg(1) = .true.
            write (output_unit, *) " Print out the physical constants used?"
            yes = read_yes_no()
            if (yes) ldbpg(2) = .true.
            ! Set options for radial modules
            write (output_unit, *) " Printout from radial modules?"
            yes = read_yes_no()
            if (yes) then
                write (output_unit, *) " Printout from RADGRD?"
                yes = read_yes_no()
                if (yes) rmcdhf_debug_flags(1) = .true.
                write (output_unit, *) " Printout from NUCPOT?"
                yes = read_yes_no()
                if (yes) rmcdhf_debug_flags(2) = .true.
                write (output_unit, *) " Printout from LODRWF?"
                yes = read_yes_no()
                if (yes) rmcdhf_debug_flags(3) = .true.
            end if
            ! Set options for angular modules
            write (output_unit, *) " Printout from angular modules?"
            yes = read_yes_no()
            if (yes) then
                write (output_unit, *) " Printout from LODCSL?"
                yes = read_yes_no()
                if (yes) ldbpa(1) = .true.
                write (output_unit, *) " Print out V coefficients?"
                yes = read_yes_no()
                if (yes) ldbpa(2) = .true.
            end if
        end if
        return
    end subroutine configure_debug_output
end module configure_debug_output_interface
