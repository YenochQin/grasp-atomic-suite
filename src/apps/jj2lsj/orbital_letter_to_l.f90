module orbital_letter_to_l_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: orbital_letter_to_l
contains

    !
    !***********************************************************************
    !                                                                      *
    integer function orbital_letter_to_l (symbol)
        !                                                                      *
        !   Modified by Gediminas Gaigalas,                   September 1997   *
        !   Modified by G. Gaigalas,                                May 2011   *
        !                                                                      *
        !***********************************************************************
        !...Switches:
        !-----------------------------------------------
        implicit none (type, external)

        ! Arguments
        character , intent(in) :: symbol

        ! Local declarations
        integer           :: letter_position
        character(len=26) :: orbital_letters

        ! Initial values
        !
        data orbital_letters/ "spdfghiklmnoqSPDFGHIKLMNOQ"/
        !-----------------------------------------------
        letter_position = index(orbital_letters,symbol)
        if (letter_position <= 13) then
            orbital_letter_to_l = letter_position - 1
        else
            orbital_letter_to_l = letter_position - 14
        end if
        return
    end function orbital_letter_to_l
end module orbital_letter_to_l_interface
