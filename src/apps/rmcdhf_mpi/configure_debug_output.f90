! Set diagnostic output options for the calculation.
module configure_debug_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_debug_output
contains

    ! This subroutine sets the arrays that control debug printout from
    ! the radial and angular modules of the GRASP92 suite.
    !
    ! Written by Farid A Parpia. Last update: 10 Dec 1992
    ! Modified bu Xinghong He. Last update: 06 Jul 1998
    subroutine configure_debug_output(debug_file_path)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit
        use debug_output_settings, only: ldbpa, ldbpg, rmcdhf_debug_flags
        use calculation_input_defaults, only: nondefault_options_flag
        use standard_io_units, only: stderr_unit
        use yes_no_input, only: read_yes_no
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        character(len=*) , intent(in) :: debug_file_path

        ! Local constants
        character(len=9), parameter :: debug_file_format = "FORMATTED"
        character(len=3), parameter :: debug_file_status = "NEW"

        ! Local declarations
        integer :: debug_path_length, specified_path_length, open_status
        character(len = len (debug_file_path)) :: file_path
        ! Initialize the diagnostic output flags.
        ! These serve as the default settings.
        ldbpa = .false.
        ldbpg = .false.

        rmcdhf_debug_flags = .false.
        if (nondefault_options_flag == 0) return

        ! Even in non-default, the user can choose not to have debug
        ! print-out

        write (stderr_unit, "(A)") "Generate debug output?  (y/n) "
        if (.not.read_yes_no()) return

        debug_path_length = len_trim(debug_file_path)

        write (stderr_unit, *) "File  ", debug_file_path(1:debug_path_length), &
            "  will be                           created as the RSCF92 DeBuG Printout&
            & File;"
        write (stderr_unit, *) "enter another file name if this is not", &
            " acceptable; null otherwise:"
        do
            read (*, "(A)") file_path

            file_path = adjustl(file_path)
            specified_path_length = len_trim(file_path)
            if (specified_path_length == 0) then
                file_path = debug_file_path
            else if (specified_path_length > debug_path_length) then
                write (stderr_unit,*) "File name too long, (> ", debug_path_length, "); redo..."
                cycle
            end if

            call open_calculation_file (angular_debug_unit, file_path, debug_file_format, debug_file_status, open_status)
            if (open_status /= 0) then
                write (stderr_unit, *) "File name not accepted; redo..."
            else
                exit
            end if
        end do
        ! Read the general output options.
        write (stderr_unit, *) "Print out the machine constants used?"
        ldbpg(1) = read_yes_no()
        write (stderr_unit, *) "Print out the physical constants used?"
        ldbpg(2) = read_yes_no()
        write (stderr_unit, *) "Printout from FNDBLK?"
        ldbpg(3) = read_yes_no()
        write (stderr_unit, *) "Print out the Hamiltonian matrix?"
        ldbpg(4) = read_yes_no()
        write (stderr_unit, *) "Print out the eigenvectors?"
        ldbpg(5) = read_yes_no()
        ! Set options for printout from radial modules
        write (stderr_unit, *) "Printout from RADGRD?"
        rmcdhf_debug_flags(1) = read_yes_no()
        write (stderr_unit, *) "Printout from NUCPOT?"
        rmcdhf_debug_flags(2) = read_yes_no()
        write (stderr_unit, *) "Printout from LODRWF?"
        rmcdhf_debug_flags(3) = read_yes_no()
        write (stderr_unit, *) "Print out I(ab) integrals?"
        rmcdhf_debug_flags(4) = read_yes_no()
        write (stderr_unit, *) "Print out Slater integrals?"
        rmcdhf_debug_flags(10) = read_yes_no()
        write (stderr_unit, *) "Make summary printout on progress", &
            " of each iteration in SOLVE?"
        rmcdhf_debug_flags(22) = read_yes_no()
        write (stderr_unit, *) "Tabulate and make printer plots", &
            " of subshell radial functions on", " each iteration in SOLVE?"
        rmcdhf_debug_flags(23) = read_yes_no()
        write (stderr_unit, *) "Tabulate and make printer plots", &
            " of subshell radial functions", " after each SCF cycle?"
        rmcdhf_debug_flags(24) = read_yes_no()
        write (stderr_unit, *) "Tabulate and make printer plots", &
            " of subshell radial functions on", " convergence?"
        rmcdhf_debug_flags(25) = read_yes_no()
        write (stderr_unit, *) "List compositions of exchange", " potentials?"
        rmcdhf_debug_flags(27) = read_yes_no()
        write (stderr_unit, *) "Tabulate and make printer plots", &
            " of exchange potentials?"
        rmcdhf_debug_flags(28) = read_yes_no()
        write (stderr_unit, *) "List compositions of direct", " potentials?"
        rmcdhf_debug_flags(29) = read_yes_no()
        write (stderr_unit, *) "Tabulate and make printer plots", &
            " of direct potentials?"
        rmcdhf_debug_flags(30) = read_yes_no()
        ! Set options for printout of angular coefficients
        write (stderr_unit, *) " Printout from LODCSL?"
        ldbpa(1) = read_yes_no()
        write (stderr_unit, *) " Print out T coefficients?"
        ldbpa(2) = read_yes_no()
        write (stderr_unit, *) " Print out V coefficients?"
        ldbpa(3) = read_yes_no()

        return
    end subroutine configure_debug_output
end module configure_debug_output_interface
