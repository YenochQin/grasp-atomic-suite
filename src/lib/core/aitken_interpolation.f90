! Interpolate a scalar table from nearby nodes using Aitken's recurrence.
module aitken_interpolation
    implicit none (type, external)

    ! Module interface
    private
    public :: interpolate_aitken
contains
    subroutine interpolate_aitken(tabulated_x, tabulated_y, table_size, target_x, interpolated_y, relative_tolerance)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: table_size
        real(real64), intent(in) :: target_x
        real(real64), intent(out) :: interpolated_y
        real(real64), intent(in) :: relative_tolerance
        real(real64), intent(in) :: tabulated_x(table_size)
        real(real64), intent(in) :: tabulated_y(table_size)

        ! Local declarations
        integer, parameter :: maximum_stencil_size = 11
        integer :: lower_bracket_index, upper_bracket_index, search_or_table_index, search_lower_bound, search_upper_bound, &
            & stencil_window_offset, stencil_row, nearest_unused_index, current_cell, diagonal_cell&
            , previous_cell, best_estimate_index, stencil_size
        real(real64) :: node_offsets(maximum_stencil_size), stencil_nodes(maximum_stencil_size), estimates(maximum_stencil_size)
        real(real64) :: aitken_table((maximum_stencil_size*(maximum_stencil_size + 1))/2)
        real(real64) :: nearest_offset, candidate_offset, best_relative_change, estimate_relative_change
        logical :: selected_nodes(2*maximum_stencil_size + 2)
        if (table_size < 1) error stop "INTERP: empty table"
        stencil_size = min(maximum_stencil_size,table_size)
        lower_bracket_index = 1
        if (target_x < tabulated_x(1)) then
            lower_bracket_index = 1
            upper_bracket_index = 1
            write (*, '("INTERP: Extrapolating, not interpolating.")')
        else if (target_x > tabulated_x(table_size)) then
            lower_bracket_index = table_size
            upper_bracket_index = table_size
            write (*, '("INTERP: Extrapolating, not interpolating.")')
        else
            search_or_table_index = 0
            do
                search_or_table_index = search_or_table_index + 1
                if (.not.(tabulated_x(search_or_table_index) < target_x)) exit
                lower_bracket_index = search_or_table_index
            end do
            upper_bracket_index = search_or_table_index
        end if
        if (target_x == tabulated_x(upper_bracket_index)) then
            interpolated_y = tabulated_y(upper_bracket_index)
            return
        end if
        search_lower_bound = max(lower_bracket_index - maximum_stencil_size,1)
        search_upper_bound = min(upper_bracket_index + maximum_stencil_size,table_size)
        stencil_window_offset = search_lower_bound - 1
        selected_nodes(search_lower_bound-stencil_window_offset:search_upper_bound-stencil_window_offset) = .false.
        ! Expand the search window as the table grows and select each node only once.
        do stencil_row = 1, stencil_size
            search_lower_bound = max(lower_bracket_index - stencil_row + 1,1)
            search_upper_bound = min(upper_bracket_index + stencil_row - 1,table_size)
            nearest_unused_index = search_lower_bound
            do while (nearest_unused_index <= search_upper_bound)
                if (.not.selected_nodes(nearest_unused_index-stencil_window_offset)) exit
                nearest_unused_index = nearest_unused_index + 1
            end do
            if (nearest_unused_index > search_upper_bound) error stop "INTERP: exhausted stencil"
            nearest_offset = tabulated_x(nearest_unused_index) - target_x
            do search_or_table_index = nearest_unused_index + 1, search_upper_bound
                if (selected_nodes(search_or_table_index-stencil_window_offset)) cycle
                candidate_offset = tabulated_x(search_or_table_index) - target_x
                if (abs(candidate_offset) < abs(nearest_offset)) then
                    nearest_offset = candidate_offset
                    nearest_unused_index = search_or_table_index
                end if
            end do
            selected_nodes(nearest_unused_index-stencil_window_offset) = .true.
            stencil_nodes(stencil_row) = tabulated_x(nearest_unused_index)
            node_offsets(stencil_row) = nearest_offset
            aitken_table(packed_triangle_index(stencil_row,1)) = tabulated_y(nearest_unused_index)
            do search_or_table_index = 2, stencil_row
                current_cell = packed_triangle_index(stencil_row,search_or_table_index)
                diagonal_cell = packed_triangle_index(search_or_table_index - 1,search_or_table_index - 1)
                previous_cell = packed_triangle_index(stencil_row,search_or_table_index - 1)
                aitken_table(current_cell) = (aitken_table(diagonal_cell)*node_offsets(stencil_row)- &
                    & aitken_table(previous_cell)*node_offsets(search_or_table_index-1))/(stencil_nodes(&
                    stencil_row)-stencil_nodes(search_or_table_index-1))
            end do
            diagonal_cell = packed_triangle_index(stencil_row,stencil_row)
            estimates(stencil_row) = aitken_table(diagonal_cell)
        end do
        if (stencil_size == 1) then
            interpolated_y = estimates(1)
            write (*,"(A)") "INTERP: One-point extrapolation; accuracy cannot be estimated."
            return
        end if
        ! Select the diagonal estimate with the smallest change relative to the previous estimate.
        best_relative_change = relative_change(estimates(2),estimates(1))

        best_estimate_index = 2
        do stencil_row = 3, stencil_size
            estimate_relative_change = relative_change(estimates(stencil_row),estimates(stencil_row-1))

            if (estimate_relative_change < best_relative_change) then
                best_relative_change = estimate_relative_change
                best_estimate_index = stencil_row
            end if

        end do
        interpolated_y = estimates(best_estimate_index)
        if (best_relative_change > relative_tolerance) write (*, '("INTERP: Accuracy of interpolation (",1p,1d10.3,&
            &") is"," below input criterion (",1d10.3,").")') best_relative_change, relative_tolerance
        return
        return
    contains

        real(real64) function relative_change(current, previous) result(change)
            ! Arguments
            real(real64), intent(in) :: current, previous
            if (current == 0.0_real64) then
                change = huge(1.0_real64)
                if (previous == 0.0_real64) change = 0.0_real64
            else
                change = abs((current-previous)/current)
            end if
        end function relative_change

        integer function packed_triangle_index (row_index, column_index)
            ! Arguments
            integer, intent(in) :: row_index
            integer, intent(in) :: column_index
            packed_triangle_index = (row_index*(row_index - 1))/2 + column_index
            return
        end function packed_triangle_index
    end subroutine interpolate_aitken
end module aitken_interpolation
