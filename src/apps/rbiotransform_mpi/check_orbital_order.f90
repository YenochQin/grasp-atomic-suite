module check_orbital_order_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: check_orbital_order
contains

    !***********************************************************************
    !                                                                      *
    subroutine check_orbital_order(n)
        !                                                                      *
        !   THIS ROUTINE DOES NOTHING!

        !   This subroutine checks the orbital order. Normal and reversed      *
        !   order can be used. If reveresed order JA and JB must be            *
        !   permuted in the subroutine ti1tv                                   *
        !                                                                      *
        !   Written by Per Jonsson                Last revision: Feb    1997   *
        !                                                                      *
        !***********************************************************************
        implicit none (type, external)

        ! Arguments
        integer  :: n
        !-----------------------------------------------

        return
    end subroutine check_orbital_order
end module check_orbital_order_interface
