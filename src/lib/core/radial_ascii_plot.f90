! Plot two tabulated radial components on the diagnostic stream using fixed-width text.
module radial_ascii_plot
    implicit none (type, external)

    ! Module interface
    private
    public :: plot_radial_orbital_ascii
contains
    subroutine plot_radial_orbital_ascii(large_component, large_component_scale, small_component, small_component_scale, &
        & last_radial_index)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use radial_grid_state, only: radial_point_count, radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: last_radial_index
        real(real64), intent(in) :: large_component_scale
        real(real64), intent(in) :: small_component_scale
        real(real64), intent(in) :: large_component(max_radial_points)
        real(real64), intent(in) :: small_component(max_radial_points)

        ! Local declarations
        integer :: plot_index, zero_baseline_row, small_component_row, large_component_row, radius_column
        real(real64) :: sixty, plot_height_scale, plot_width_scale, maximum_amplitude, minimum_amplitude, &
            & scaled_large_value, scaled_small_value, maximum_radius, amplitude_to_row, &
            index_to_column, small_component_to_row, large_component_to_row
        logical :: initialize_plot_templates
        character :: bordered_row*132, dashed_row*132
        character  :: plot_rows(last_radial_index)*132
        character :: large_component_marker, small_component_marker

        ! Initial values
        data sixty/ 60.0_real64/
        data plot_height_scale/ 54.0_real64/
        data plot_width_scale/ 131.0_real64/
        data initialize_plot_templates/ .true./

        if (initialize_plot_templates) then
            bordered_row(1:1) = "|"
            dashed_row(1:1) = "-"
            do plot_index = 2, 131
                bordered_row(plot_index:plot_index) = " "
                dashed_row(plot_index:plot_index) = "-"
            end do
            bordered_row(132:132) = "|"
            dashed_row(132:132) = "-"
            initialize_plot_templates = .false.
        end if
        if (small_component_scale == 0.0_real64) then
            large_component_marker = "X"
            small_component_marker = "-"
        else
            large_component_marker = "P"
            small_component_marker = "Q"
        end if
        maximum_amplitude = max(large_component_scale*large_component(1),small_component_scale*small_component(1))
        minimum_amplitude = min(large_component_scale*large_component(1),small_component_scale*small_component(1))
        do plot_index = 2, last_radial_index
            scaled_large_value = large_component_scale*large_component(plot_index)
            scaled_small_value = small_component_scale*small_component(plot_index)
            maximum_amplitude = max(max(maximum_amplitude,scaled_large_value),scaled_small_value)
            minimum_amplitude = min(min(minimum_amplitude,scaled_large_value),scaled_small_value)
        end do
        maximum_radius = radial_grid(last_radial_index)
        if (maximum_amplitude==0.0_real64 .and. minimum_amplitude==0.0_real64) then
            call ensure_debug_output()
            write (angular_debug_unit, '(" No plot: function is identically zero")')
            return
        else
            amplitude_to_row = plot_height_scale/abs(maximum_amplitude - minimum_amplitude)
        end if
        index_to_column = plot_width_scale/dble(radial_point_count)
        if (maximum_amplitude>0.0_real64 .and. minimum_amplitude<0.0_real64) then
            zero_baseline_row = 4 + amplitude_to_row*abs(minimum_amplitude)
        else
            zero_baseline_row = 1
        end if
        do plot_index = 1, 60
            if (plot_index==1 .or. plot_index==zero_baseline_row .or. plot_index==60) then
                plot_rows(plot_index) = dashed_row
            else
                plot_rows(plot_index) = bordered_row
            end if
        end do
        small_component_to_row = small_component_scale*amplitude_to_row
        large_component_to_row = large_component_scale*amplitude_to_row
        do plot_index = 1, last_radial_index
            small_component_row = nint(small_component_to_row*small_component(plot_index)) + zero_baseline_row
            large_component_row = nint(large_component_to_row*large_component(plot_index)) + zero_baseline_row
            radius_column = max(1,nint(dble(plot_index)*index_to_column))
            plot_rows(small_component_row)(radius_column:radius_column) = small_component_marker
            plot_rows(large_component_row)(radius_column:radius_column) = large_component_marker
        end do
        call ensure_debug_output()
        write (angular_debug_unit, '("1")')
        do plot_index = 60, 1, -1
            call ensure_debug_output()
            write (angular_debug_unit, '(1x,1a132)') plot_rows(plot_index)
        end do
        call ensure_debug_output()
        write (angular_debug_unit, '(/,50x,1p," r (max) = ",1d10.3,"Bohr radii,"," Maximum of functions = ",1d10.3)') &
            maximum_radius, maximum_amplitude
        return
        return
    end subroutine plot_radial_orbital_ascii
end module radial_ascii_plot
