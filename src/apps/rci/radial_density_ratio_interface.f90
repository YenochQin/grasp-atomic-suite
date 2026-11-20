module radial_density_ratio_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function radial_density_ratio (large_component_values, small_component_values, last_radial_index, &
            & principal_n, kappa, &
            & effective_nuclear_charge)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,   only: max_radial_points
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: large_component_values(max_radial_points)
            real(real64), intent(in) :: small_component_values(max_radial_points)
            integer, intent(in) :: last_radial_index
            integer :: principal_n
            integer :: kappa
            real(real64) :: effective_nuclear_charge
        end function radial_density_ratio
    end interface
end module radial_density_ratio_interface
