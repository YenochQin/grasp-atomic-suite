module specific_mass_shift_radial_integral_3_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        real(real64) function specific_mass_shift_radial_integral_3 (first_orbital, second_orbital)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: first_orbital
            integer, intent(in) :: second_orbital
        end function specific_mass_shift_radial_integral_3
    end interface
end module specific_mass_shift_radial_integral_3_interface
