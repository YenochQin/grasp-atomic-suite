module normal_mass_shift_radial_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine normal_mass_shift_radial_integral (first_orbital, second_orbital, relativistic_nms_integral,leading_nms_integral)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(out) :: relativistic_nms_integral, leading_nms_integral
            integer, intent(in) :: first_orbital
            integer, intent(in) :: second_orbital
        end subroutine normal_mass_shift_radial_integral
    end interface
end module normal_mass_shift_radial_integral_interface
