! Integrate a radial integrand using an origin power-law fit and composite five-point quadrature.
module radial_quadrature
    implicit none (type, external)

    ! Module interface
    private
    public :: integrate_radial_grid
contains
    subroutine integrate_radial_grid(integral)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_grid, radial_grid_derivative
        use radial_quadrature_weights, only: endpoint_weight, odd_interior_weight, middle_weight, shared_endpoint_weight
        use radial_quadrature_workspace, only: radial_workspace_a, integration_last_index
        implicit none (type, external)

        ! Arguments
        real(real64) , intent(out) :: integral

        ! Local declarations
        integer :: penultimate_grid_index, radial_index, next_radial_index
        real(real64) :: weighted_value, next_weighted_value, weighted_value_ratio, next_integrand, integrand, &
            & integrand_ratio, next_radius, radius, origin_power
        ! Workspace A contains f(r) times dr/dt. Find a usable pair for the small-radius power fit.
        penultimate_grid_index = integration_last_index - 1
        do radial_index = 2, penultimate_grid_index
            weighted_value = radial_workspace_a(radial_index)
            if (abs(weighted_value) <= 0.0_real64) cycle
            next_radial_index = radial_index + 1
            next_weighted_value = radial_workspace_a(next_radial_index)
            weighted_value_ratio = next_weighted_value/weighted_value
            if (weighted_value_ratio <= 0.0_real64) cycle
            next_integrand = next_weighted_value/radial_grid_derivative(next_radial_index)
            integrand = weighted_value/radial_grid_derivative(radial_index)
            integrand_ratio = next_integrand/integrand
            next_radius = radial_grid(next_radial_index)
            radius = radial_grid(radial_index)
            origin_power = log(integrand_ratio)/log(next_radius/radius)
            integrand = radius*integrand
            integral = integrand/(origin_power + 1.0_real64)
            ! Pad the tail in the shared workspace for the last five-point panel.
            radial_workspace_a(integration_last_index+1:3+integration_last_index) = 0.0_real64
            integral = integral + endpoint_weight*weighted_value
            integral = integral + sum(odd_interior_weight*(radial_workspace_a(next_radial_index:integration_last_index:4) + &
                radial_workspace_a(next_radial_index+2:integration_last_index+2:4)) + &
                middle_weight*radial_workspace_a(next_radial_index+1:integration_last_index+1:4) + &
                shared_endpoint_weight*radial_workspace_a(next_radial_index+3:integration_last_index+3:4))
            ! A final point on a panel boundary received the weight of a shared endpoint.
            if (mod(integration_last_index - radial_index, 4) == 0) &
                integral = integral - endpoint_weight*radial_workspace_a(integration_last_index)
            return
        end do
        integral = 0.0_real64
    end subroutine integrate_radial_grid
end module radial_quadrature
