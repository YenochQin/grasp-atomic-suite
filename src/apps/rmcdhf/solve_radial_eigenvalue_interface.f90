module solve_radial_eigenvalue_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function solve_radial_eigenvalue (orbital_index)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
        end function solve_radial_eigenvalue
    end interface
end module solve_radial_eigenvalue_interface
