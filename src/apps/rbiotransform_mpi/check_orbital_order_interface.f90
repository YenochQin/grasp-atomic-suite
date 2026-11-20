module check_orbital_order_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine check_orbital_order (n)
            implicit none (type, external)

            ! Arguments
            integer :: n
        end subroutine check_orbital_order
    end interface
end module check_orbital_order_interface
