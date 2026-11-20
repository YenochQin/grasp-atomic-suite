module calculate_vacuum_polarization_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine calculate_vacuum_polarization_integral (first_orbital,second_orbital,radial_integral_value)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer , intent(inout) :: first_orbital, second_orbital
            real(real64), intent(out) :: radial_integral_value
        end subroutine calculate_vacuum_polarization_integral
    end interface
end module calculate_vacuum_polarization_integral_interface
