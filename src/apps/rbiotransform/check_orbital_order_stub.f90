! Retain the historical orbital-order entry point as a stub that returns immediately.
module check_orbital_order_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: check_orbital_order
contains

    ! Written by Per Jonsson. Last revision: Feb. 1997
    subroutine check_orbital_order(n)

        implicit none (type, external)

        ! Arguments
        integer  :: n

        return
    end subroutine check_orbital_order
end module check_orbital_order_interface
