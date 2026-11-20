module field_shift_radial_integral_vector_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine field_shift_radial_integral_vector (first_orbital, second_orbital, dint1vec, nuclear_grid_point_count)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,    only: max_orbitals
            use radial_grid_state,           only: radial_point_count
            implicit none (type, external)

            ! Arguments
            real(real64), intent(out) :: dint1vec(max_orbitals, max_orbitals, radial_point_count)
            integer, intent(in) :: first_orbital
            integer, intent(in) :: second_orbital
            integer, intent(in) :: nuclear_grid_point_count
        end subroutine field_shift_radial_integral_vector
    end interface
end module field_shift_radial_integral_vector_interface
