module exchange_angular_coefficient_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function exchange_angular_coefficient (multipole_rank, csf_index, orbital_a, orbital_b)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer :: multipole_rank
            integer :: csf_index
            integer, intent(in) :: orbital_a
            integer, intent(in) :: orbital_b
        end function exchange_angular_coefficient
    end interface
end module exchange_angular_coefficient_interface
