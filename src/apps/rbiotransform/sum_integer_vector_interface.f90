module sum_integer_vector_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function sum_integer_vector (ivec, nelmnt)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ivec(*)
            integer, intent(in) :: nelmnt
        end function sum_integer_vector
    end interface
end module sum_integer_vector_interface
