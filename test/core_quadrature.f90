program core_quadrature
    implicit none (type, external)

    call setup
    call runtest

contains

    subroutine runtest
        use grasp_parameters, only: radial_workspace_capacity
        use radial_grid_state, only: transformed_grid_step, asymptotic_grid_step, radial_point_count, &
            radial_workspace_capacity, radial_grid, radial_grid_scale, radial_grid_derivative, radial_derivative_over_radius
        use radial_quadrature_workspace, only: radial_workspace_capacity, integration_last_index, radial_workspace_a
        use radial_quadrature, only: integrate_radial_grid
        use, intrinsic :: iso_fortran_env, only: real64


        integer :: i
        real(real64) :: result

        ! Integrate exp(-2x) from 0 to infinity. Expected result is 1/2.
        do i = 1, radial_workspace_capacity
            radial_workspace_a(i) = exp(-2*radial_grid(i)) * radial_grid_derivative(i)
        end do
        call integrate_radial_grid(result)

        print *, "RESULT from QUAD", result
        print *, "expected        ", 0.5_real64
        print *, "difference      ", abs(result - 0.5_real64)

        if(abs(result - 0.5_real64) > 1e-12_real64) then
            stop 1
        end if
    end subroutine runtest

    subroutine setup
        use grasp_parameters, only: max_radial_points
        use radial_grid_state, only: transformed_grid_step, asymptotic_grid_step, radial_point_count, &
            radial_workspace_capacity, radial_grid, radial_grid_scale, radial_grid_derivative, radial_derivative_over_radius
        use radial_quadrature_workspace, only: radial_workspace_capacity, integration_last_index, radial_workspace_a
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        use radial_grid_setup, only: initialize_radial_grid
        use, intrinsic :: iso_fortran_env, only: real64

        transformed_grid_step = 5.0e-2_real64
        radial_grid_scale = 2.0e-6_real64
        asymptotic_grid_step = 0.0_real64
        radial_point_count = max_radial_points
        integration_last_index = max_radial_points

        call initialize_radial_stencil_weights
        call initialize_radial_grid
    end subroutine setup

end program core_quadrature
