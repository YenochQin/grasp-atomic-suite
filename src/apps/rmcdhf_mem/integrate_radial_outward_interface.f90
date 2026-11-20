module integrate_radial_outward_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine integrate_radial_outward (orbital_index, matching_grid_index, large_component_values, small_component_values)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,    only: max_radial_points
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            integer, intent(in) :: matching_grid_index
            real(real64), intent(inout) :: large_component_values(max_radial_points)
            real(real64), intent(inout) :: small_component_values(max_radial_points)
        end subroutine integrate_radial_outward
    end interface
end module integrate_radial_outward_interface
