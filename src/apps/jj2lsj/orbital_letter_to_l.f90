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
    integer           :: locate
    character(len=26) :: set

    ! Initial values
    !
    data set/ "spdfghiklmnoqSPDFGHIKLMNOQ"/
    !-----------------------------------------------
    locate = index(set,symbol)
    if (locate <= 13) then
        orbital_letter_to_l = locate - 1
    else
        orbital_letter_to_l = locate - 14
    end if
    return
end function orbital_letter_to_l
