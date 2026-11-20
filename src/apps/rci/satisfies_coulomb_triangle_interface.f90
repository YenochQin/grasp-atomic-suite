module satisfies_coulomb_triangle_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        logical function satisfies_coulomb_triangle (la, k, lb)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: la
            integer, intent(in) :: k
            integer, intent(in) :: lb
        end function satisfies_coulomb_triangle
    end interface
end module satisfies_coulomb_triangle_interface
