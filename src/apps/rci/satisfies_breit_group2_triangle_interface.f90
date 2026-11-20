module satisfies_breit_group2_triangle_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        logical function satisfies_breit_group2_triangle (orbital_a,orbital_b,orbital_c,orbital_d,l)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_a,orbital_b,orbital_c,orbital_d,l
        end function satisfies_breit_group2_triangle
    end interface
end module satisfies_breit_group2_triangle_interface
