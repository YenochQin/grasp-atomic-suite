module find_least_consistent_orbital_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine find_least_consistent_orbital (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(out) :: orbital_index
        end subroutine find_least_consistent_orbital
    end interface
end module find_least_consistent_orbital_interface
