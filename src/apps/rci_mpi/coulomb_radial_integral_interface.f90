module coulomb_radial_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function coulomb_radial_integral (recoupling_coefficient, orbital_a, orbital_c, rbd, orbital_b, orbital_d, &
            & multipole_rank, iw)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,   only: max_radial_points
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: recoupling_coefficient(max_radial_points)
            integer, intent(in) :: orbital_a
            integer, intent(in) :: orbital_c
            real(real64), intent(in) :: rbd(max_radial_points)
            integer, intent(in) :: orbital_b
            integer, intent(in) :: orbital_d
            integer, intent(in) :: multipole_rank
            integer, intent(in) :: iw
        end function coulomb_radial_integral
    end interface
end module coulomb_radial_integral_interface
