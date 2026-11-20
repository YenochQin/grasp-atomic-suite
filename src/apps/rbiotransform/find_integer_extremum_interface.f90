module find_integer_extremum_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function find_integer_extremum (ivec, nel, imxmn)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nel
            integer, intent(in) :: ivec(nel)
            integer, intent(in) :: imxmn
            !...This routine performs I/O.
        end function find_integer_extremum
    end interface
end module find_integer_extremum_interface
