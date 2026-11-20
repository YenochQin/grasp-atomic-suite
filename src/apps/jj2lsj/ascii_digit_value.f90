!***********************************************************************
!                                                                      *
integer function ascii_digit_value (cst)
    !                                                                      *
    !                                                                      *
    !                                                                      *
    !***********************************************************************
    !-----------------------------------------------
    implicit none (type, external)

    ! Arguments
    character , intent(in) :: cst

    ! Local declarations
    integer :: i
    character :: cdgt(0:9)

    ! Initial values
    !-----------------------------------------------
    !
    data cdgt/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"/
    !
    do i = 0, 9
        if (cst /= cdgt(i)) cycle
        ascii_digit_value = i
        exit
    end do
    !
    return
end function ascii_digit_value
