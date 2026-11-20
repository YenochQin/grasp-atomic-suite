! Transfer both components of a radial orbital to the current grid and renormalize them.
module radial_orbital_interpolation
    implicit none (type, external)

    ! Module interface
    private
    public :: interpolate_radial_orbital
contains
    subroutine interpolate_radial_orbital(source_large_component, source_small_component, source_grid_size, source_radial_grid, &
        & orbital_index, original_norm_squared)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state,           only: working_tolerance
        use radial_grid_state, only: radial_point_count, radial_grid
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: source_grid_size
        integer, intent(in) :: orbital_index
        real(real64), intent(out) :: original_norm_squared
        real(real64), intent(in) :: source_large_component(*), source_small_component(*), source_radial_grid(*)

        ! Local declarations
        integer, parameter :: maximum_stencil_size = 13
        integer :: radial_index, last_target_index, lower_bracket_index, unconverged_point_count, stencil_row, &
            & search_or_table_index, upper_bracket_index, search_lower_bound, search_upper_bound, nearest_unused_index, &
            current_cell, diagonal_cell, previous_cell, first_tail_index
        real(real64) :: source_maximum_radius
        real(real64) :: stencil_nodes(maximum_stencil_size), node_offsets(maximum_stencil_size)
        real(real64) :: large_component_table((maximum_stencil_size*(maximum_stencil_size + 1))/2), &
            & small_component_table((maximum_stencil_size*(maximum_stencil_size + 1))/2)
        real(real64) :: target_maximum_radius, target_radius, previous_large_estimate, previous_small_estimate, &
            & nearest_offset, candidate_offset, previous_node_offset, current_node_offset, &
            inverse_node_separation, large_estimate, small_estimate, large_relative_change, small_relative_change, &
            & normalization_factor
        logical :: found_candidate
        logical :: selected_nodes(max_radial_points)
        source_maximum_radius = source_radial_grid(source_grid_size)
        target_maximum_radius = radial_grid(radial_point_count)
        radial_large_components(1,orbital_index) = 0.0_real64
        radial_small_components(1,orbital_index) = 0.0_real64
        ! The current grid must reach the end of the input grid.
        if (source_maximum_radius > target_maximum_radius) then
            write (*, '(/,"INTRPQ: Grid of insufficient extent:"/," Present grid has R(N) = ",1p,1d19.12,&
                &" Bohr radii"/,"          Require R(N) = ",1d19.12," Bohr radii")') &
                target_maximum_radius, source_maximum_radius
            print*, "N =",radial_point_count,"MA =",source_grid_size
            stop
        end if
        radial_index = radial_point_count
        do
            radial_index = radial_index - 1
            if (radial_grid(radial_index) <= source_maximum_radius) exit
        end do
        last_target_index = radial_index
        orbital_last_grid_index(orbital_index) = last_target_index
        lower_bracket_index = 0
        unconverged_point_count = 0
        do radial_index = 2, last_target_index
            target_radius = radial_grid(radial_index)
            stencil_row = 0
            previous_large_estimate = 0.0_real64
            previous_small_estimate = 0.0_real64
            do
                search_or_table_index = lower_bracket_index + 1
                if (.not.(source_radial_grid(search_or_table_index) < target_radius)) exit
                lower_bracket_index = search_or_table_index
            end do
            upper_bracket_index = search_or_table_index
            search_lower_bound = max(lower_bracket_index - maximum_stencil_size,1)
            search_upper_bound = min(upper_bracket_index + maximum_stencil_size,source_grid_size)
            selected_nodes(search_lower_bound:search_upper_bound) = .false.
            do
                stencil_row = stencil_row + 1
                search_lower_bound = max(lower_bracket_index - stencil_row + 1,1)
                search_upper_bound = min(upper_bracket_index + stencil_row - 1,source_grid_size)
                found_candidate = .false.
                do search_or_table_index = search_lower_bound, search_upper_bound
                    if (selected_nodes(search_or_table_index)) cycle
                    if (.not.found_candidate) then
                        nearest_offset = source_radial_grid(search_or_table_index) - target_radius
                        nearest_unused_index = search_or_table_index
                        found_candidate = .true.
                    else
                        candidate_offset = source_radial_grid(search_or_table_index) - target_radius
                        if (abs(candidate_offset) < abs(nearest_offset)) then
                            nearest_offset = candidate_offset
                            nearest_unused_index = search_or_table_index
                        end if
                    end if
                end do
                ! Add the same nearby node to both Aitken tables. Each diagonal gives the estimate
                ! at the target radius using all nodes selected so far.
                selected_nodes(nearest_unused_index) = .true.
                stencil_nodes(stencil_row) = source_radial_grid(nearest_unused_index)
                node_offsets(stencil_row) = nearest_offset
                do search_or_table_index = 1, stencil_row
                    current_cell = packed_triangle_index(stencil_row,search_or_table_index)
                    if (search_or_table_index == 1) then
                        large_component_table(current_cell) = source_large_component(nearest_unused_index)
                        small_component_table(current_cell) = source_small_component(nearest_unused_index)
                    else
                        diagonal_cell = packed_triangle_index(search_or_table_index - 1,search_or_table_index - 1)
                        previous_cell = packed_triangle_index(stencil_row,search_or_table_index - 1)
                        previous_node_offset = node_offsets(search_or_table_index-1)
                        current_node_offset = node_offsets(stencil_row)
                        inverse_node_separation = 1.0_real64/(stencil_nodes(stencil_row)-stencil_nodes(search_or_table_index-1))
                        large_component_table(current_cell) = &
                            & (large_component_table(diagonal_cell)*current_node_offset-large_component_table(previous_cell)* &
                            & previous_node_offset)*&
                            inverse_node_separation
                        small_component_table(current_cell) = &
                            & (small_component_table(diagonal_cell)*current_node_offset-small_component_table(previous_cell)* &
                            & previous_node_offset)*&
                            inverse_node_separation
                    end if
                end do
                diagonal_cell = packed_triangle_index(stencil_row,stencil_row)
                large_estimate = large_component_table(diagonal_cell)
                small_estimate = small_component_table(diagonal_cell)
                if (large_estimate==0.0_real64 .or. small_estimate==0.0_real64) then
                    if (stencil_row < maximum_stencil_size) cycle
                else
                    large_relative_change = abs((large_estimate - previous_large_estimate)/large_estimate)
                    small_relative_change = abs((small_estimate - previous_small_estimate)/small_estimate)
                    if (.not.(small_relative_change < working_tolerance .and. large_relative_change < working_tolerance)) then
                        previous_large_estimate = large_estimate
                        previous_small_estimate = small_estimate
                        if (stencil_row < maximum_stencil_size) cycle
                        unconverged_point_count = unconverged_point_count + 1
                    end if
                end if
                radial_large_components(radial_index,orbital_index) = large_estimate
                radial_small_components(radial_index,orbital_index) = small_estimate
                exit
            end do
        end do
        first_tail_index = last_target_index + 1
        radial_large_components(first_tail_index:radial_point_count,orbital_index) = 0.0_real64
        radial_small_components(first_tail_index:radial_point_count,orbital_index) = 0.0_real64
        if (rmcdhf_debug_flags(3) .and. unconverged_point_count>0) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,&
                &"INTRPQ: Interpolation procedure not converged to",1p,1d19.12," for ",1i3," of ",1i3,&
                &" tabulation points")') working_tolerance, unconverged_point_count, &
                & last_target_index
        end if
        ! Return the norm squared before scaling both interpolated components.
        original_norm_squared = integrate_weighted_radial_overlap(orbital_index,orbital_index,0)
        normalization_factor = 1.0_real64/sqrt(original_norm_squared)
        radial_large_components(:last_target_index,orbital_index) = radial_large_components(:last_target_index, &
            & orbital_index)*normalization_factor
        radial_small_components(:last_target_index,orbital_index) = radial_small_components(:last_target_index, &
            & orbital_index)*normalization_factor
        return
        return
    contains

        integer function packed_triangle_index (row_index, column_index)
            ! Arguments
            integer, intent(in) :: row_index
            integer, intent(in) :: column_index
            packed_triangle_index = (row_index*(row_index - 1))/2 + column_index
            return
        end function packed_triangle_index
    end subroutine interpolate_radial_orbital
end module radial_orbital_interpolation
