module solve_homogeneous_radial_equation_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine solve_homogeneous_radial_equation (orbital_index, solve_failed)
            ! Module dependencies
            use grasp_parameters, only: max_radial_points, radial_workspace_capacity, max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            logical, intent(out) :: solve_failed
        end subroutine solve_homogeneous_radial_equation
    end interface
end module solve_homogeneous_radial_equation_interface
