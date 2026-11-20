module integrate_homogeneous_radial_segment_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine integrate_homogeneous_radial_segment (orbital_index, first_radial_index, last_radial_index, &
            & large_component_values, small_component_values)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters, only: max_radial_points, radial_workspace_capacity, max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            integer, intent(in) :: first_radial_index
            integer, intent(in) :: last_radial_index
            real(real64), intent(inout) :: large_component_values(max_radial_points)
            real(real64), intent(inout) :: small_component_values(max_radial_points)
        end subroutine integrate_homogeneous_radial_segment
    end interface
end module integrate_homogeneous_radial_segment_interface
