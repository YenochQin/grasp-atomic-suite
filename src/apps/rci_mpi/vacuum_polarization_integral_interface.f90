module vacuum_polarization_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        real(real64)  function vacuum_polarization_integral (first_orbital,second_orbital)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer , intent(in) :: first_orbital, second_orbital
        end function vacuum_polarization_integral
    end interface
end module vacuum_polarization_integral_interface
