! Read a case-insensitive yes/no response, prompting again for other characters.
module yes_no_input
    implicit none (type, external)

    ! Module interface
    private
    public :: getyn
contains
    logical function getyn ()
        ! Module dependencies
        use standard_io_units,     only: stdin_unit, stderr_unit
        implicit none (type, external)

        ! Local declarations
        character (len = 1) :: response_character
        do
            read (stdin_unit, "(A)") response_character
            select case (response_character)
            case ("y", "Y")
                getyn = .true.
                return
            case ("n", "N")
                getyn = .false.
                return
            case default
                write (stderr_unit, *) "Expecting <y><cr> or <n><cr> ..."
            end select
        end do
    end function getyn
end module yes_no_input
