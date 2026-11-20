! Evaluate the inner and outer Hartree radial integrals and publish the Z and Y functions.
module hartree_yz_functions
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_hartree_yz_functions
contains
    subroutine evaluate_hartree_yz_functions(tensor_rank, first_orbital, second_orbital)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: radial_workspace_capacity
        use radial_difference_coefficients,           only: five_point_weights
        use atomic_calculation_state,           only: working_tolerance
        use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative
        use radial_quadrature_weights, only: endpoint_weight, odd_interior_weight, middle_weight
        use radial_quadrature_workspace, only: radial_workspace_a, radial_workspace_b, integration_last_index
        use radial_wavefunctions,          only: orbital_last_grid_index, radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: tensor_rank
        integer, intent(in) :: first_orbital
        integer, intent(in) :: second_orbital

        ! Local declarations
        integer, parameter :: maximum_tensor_rank = 20
        integer, parameter :: rank_cache_size = maximum_tensor_rank + 1
        integer :: radial_index, first_padding_index, last_quadrature_index, first_asymptotic_index, stencil_or_tail_index, &
            & padded_grid_end
        real(real64) :: weighted_orbital_product(radial_workspace_capacity)
        real(real64) :: radius_to_rank(radial_workspace_capacity, maximum_tensor_rank), &
            & inverse_radius_to_rank(radial_workspace_capacity, maximum_tensor_rank), &
            radius_to_rank_plus_one(radial_workspace_capacity, maximum_tensor_rank), &
            & inverse_radius_to_rank_plus_one(radial_workspace_capacity, maximum_tensor_rank)
        real(real64) :: inverse_radius(radial_workspace_capacity), hartree_outer_integral(radial_workspace_capacity), &
            & weighted_integrand(radial_workspace_capacity)
        real(real64) :: stencil_integral, interval_integral, inner_integral_limit, tail_correction
        logical :: rank_cached(0:maximum_tensor_rank)

        save rank_cached, radius_to_rank, inverse_radius_to_rank, radius_to_rank_plus_one, inverse_radius_to_rank_plus_one, &
            & inverse_radius

        ! Initial values
        data rank_cached/ rank_cache_size*.false./
        ! Cache the radial powers for each tensor rank on the initialized grid.
        ! The cache key contains only the rank, so the radial grid must remain fixed across calls.
        if (.not.rank_cached(tensor_rank)) then
            rank_cached(tensor_rank) = .true.
            if (tensor_rank > maximum_tensor_rank) then
                write (output_unit, *) " increase maximum_tensor_rank in evaluate_hartree_yz_functions "
                stop
            end if
            if (tensor_rank > 0) then
                radius_to_rank(2:radial_point_count,tensor_rank) = radial_grid(2:radial_point_count)**tensor_rank
                inverse_radius_to_rank(2:radial_point_count,tensor_rank) = 1.0_real64/radius_to_rank(2:radial_point_count, &
                    & tensor_rank)
                radius_to_rank_plus_one(2:radial_point_count,tensor_rank) = radius_to_rank(2:radial_point_count, &
                    & tensor_rank)*radial_grid(2:radial_point_count)
                inverse_radius_to_rank_plus_one(2:radial_point_count, &
                    & tensor_rank) = 1.0_real64/radius_to_rank_plus_one(2:radial_point_count,tensor_rank)
                inverse_radius(2:radial_point_count) = 1.0_real64/radial_grid(2:radial_point_count)
            else
                inverse_radius(2:radial_point_count) = 1.0_real64/radial_grid(2:radial_point_count)
            end if
        end if
        integration_last_index = min(orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        first_padding_index = integration_last_index + 1
        last_quadrature_index = integration_last_index + 3
        first_asymptotic_index = integration_last_index + 4
        do radial_index = 2, integration_last_index
            weighted_orbital_product(radial_index) = &
                & radial_grid_derivative(radial_index)*(radial_large_components(radial_index, &
                & first_orbital)*radial_large_components(radial_index,second_orbital) + &
                & radial_small_components(radial_index,first_orbital)*radial_small_components(radial_index,second_orbital))
        end do
        weighted_integrand(1) = 0.0_real64
        if (tensor_rank == 0) then
            weighted_integrand(2:integration_last_index) = weighted_orbital_product(2:integration_last_index)
        else
            weighted_integrand(2:integration_last_index) = radius_to_rank(2:integration_last_index, &
                & tensor_rank)*weighted_orbital_product(2:integration_last_index)
        end if
        weighted_integrand(first_padding_index:first_asymptotic_index) = 0.0_real64
        ! Accumulate the inner integral from the origin before applying the radial power. Store Z^k in workspace A.
        radial_workspace_a(1) = 0.0_real64
        do radial_index = 2, 4
            stencil_integral = 0.0_real64
            do stencil_or_tail_index = 2, 5
                stencil_integral = stencil_integral + five_point_weights(stencil_or_tail_index, &
                    & radial_index)*weighted_integrand(stencil_or_tail_index)
            end do
            radial_workspace_a(radial_index) = stencil_integral
        end do
        do radial_index = 5, last_quadrature_index
            interval_integral = endpoint_weight*(weighted_integrand(radial_index-4)+weighted_integrand(radial_index)) + &
                & odd_interior_weight*(weighted_integrand(radial_index-3)+weighted_integrand(radial_index-1)) + &
                & middle_weight*weighted_integrand&
                (radial_index-2)
            radial_workspace_a(radial_index) = radial_workspace_a(radial_index-4) + interval_integral
        end do
        ! For k=0, impose the orthonormal overlap limit on each of the four integration chains.
        if (tensor_rank == 0) then
            if (first_orbital == second_orbital) then
                inner_integral_limit = 1.0_real64
            else
                inner_integral_limit = 0.0_real64
            end if
            do stencil_or_tail_index = last_quadrature_index, integration_last_index, -1
                tail_correction = radial_workspace_a(stencil_or_tail_index) - inner_integral_limit
                if (abs(tail_correction) <= working_tolerance) cycle
                radial_workspace_a(stencil_or_tail_index:stencil_or_tail_index-((-2-stencil_or_tail_index)/(-4)-1)*4:(-4)) = &
                    & radial_workspace_a(stencil_or_tail_index:stencil_or_tail_index-((-2-stencil_or_tail_index)/(-4)-1)*4:&
                    (-4)) - tail_correction
            end do
        else
            inner_integral_limit = radial_workspace_a(last_quadrature_index)
        end if
        if (tensor_rank == 0) then
            radial_workspace_a(first_asymptotic_index:radial_point_count) = inner_integral_limit
        else
            radial_workspace_a(2:last_quadrature_index) = &
                & radial_workspace_a(2:last_quadrature_index)*inverse_radius_to_rank(2:last_quadrature_index,tensor_rank)
            radial_workspace_a(first_asymptotic_index:radial_point_count) = &
                & inner_integral_limit*inverse_radius_to_rank(first_asymptotic_index:radial_point_count,tensor_rank)
        end if
        ! Start the outer integral beyond the common orbital support and accumulate inwards.
        padded_grid_end = radial_point_count + 4
        hartree_outer_integral(padded_grid_end:first_padding_index:(-1)) = 0.0_real64
        weighted_integrand(1) = 0.0_real64
        if (tensor_rank == 0) then
            weighted_integrand(2:integration_last_index) = &
                & weighted_orbital_product(2:integration_last_index)*inverse_radius(2:integration_last_index)
        else
            weighted_integrand(2:integration_last_index) = &
                & weighted_orbital_product(2:integration_last_index)* &
                    inverse_radius_to_rank_plus_one(2:integration_last_index,tensor_rank)
        end if
        do radial_index = integration_last_index, 2, -1
            hartree_outer_integral(radial_index) = hartree_outer_integral(radial_index+4) + &
                & endpoint_weight*(weighted_integrand(radial_index)+weighted_integrand(radial_index+4)) + &
                & odd_interior_weight*(weighted_integrand(radial_index+1)+weighted_integrand(radial_index+&
                3)) + middle_weight*weighted_integrand(radial_index+2)
        end do
        hartree_outer_integral(1) = 0.0_real64
        if (tensor_rank == 0) then
            hartree_outer_integral(2:integration_last_index) = &
                & hartree_outer_integral(2:integration_last_index)*radial_grid(2:integration_last_index)
        else
            hartree_outer_integral(2:integration_last_index) = &
                & hartree_outer_integral(2:integration_last_index)*radius_to_rank_plus_one(2:integration_last_index,tensor_rank)
        end if
        ! Store Y^k in workspace B, adding the outer integral to Z^k.
        radial_workspace_b(1) = 0.0_real64
        radial_workspace_b(2:radial_point_count) = radial_workspace_a(2:radial_point_count) + &
            & hartree_outer_integral(2:radial_point_count)
        return
    end subroutine evaluate_hartree_yz_functions
end module hartree_yz_functions
