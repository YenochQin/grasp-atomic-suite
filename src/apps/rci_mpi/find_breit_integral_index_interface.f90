module find_breit_integral_index_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function find_breit_integral_index (itype, i)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: itype
            integer :: i
        end function find_breit_integral_index
    end interface
end module find_breit_integral_index_interface
