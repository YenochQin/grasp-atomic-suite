! Read a case-insensitive yes/no response, prompting again for other characters.
module yes_no_input
    implicit none (type, external)

    ! Module interface
    private
    public :: read_yes_no
contains
    logical function read_yes_no ()
        ! Module dependencies
        use standard_io_units,     only: stdin_unit, stderr_unit
        implicit none (type, external)

        ! Local declarations
        character (len = 1) :: response_character
        do
            read (stdin_unit, "(A)") response_character
            select case (response_character)
            case ("y", "Y")
                read_yes_no = .true.
                return
            case ("n", "N")
                read_yes_no = .false.
                return
            case default
                write (stderr_unit, *) "Expecting <y><cr> or <n><cr> ..."
            end select
        end do
    end function read_yes_no
end module yes_no_input
