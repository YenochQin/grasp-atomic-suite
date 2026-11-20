! Count sign changes between significant radial extrema and return the sign of the first retained extremum.
module radial_node_count
    implicit none (type, external)

    ! Module interface
    private
    public :: count_radial_nodes
contains
    subroutine count_radial_nodes(radial_values, last_radial_index, node_count, first_extremum_sign)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit, ensure_calculation_log
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_radial_points
        use radial_node_count_settings, only: count_context, count_trace_rank, relative_extremum_cutoff
        use atomic_calculation_state, only: working_tolerance
        use radial_grid_state, only: radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: last_radial_index
        integer, intent(out) :: node_count
        real(real64), intent(out) :: first_extremum_sign
        real(real64), intent(in) :: radial_values(max_radial_points)

        ! Local declarations
        integer :: extremum_grid_indices(max_radial_points)
        integer :: extremum_count, sample_index, extremum_index, remaining_extrema
        integer :: trace_status
        character(len=8) :: trace_setting
        logical, save :: trace_enabled = .false., trace_initialized = .false.
        real(real64) :: segment_max_amplitude, peak_amplitude, sample_amplitude, absolute_sign_sum, extremum_cutoff, &
            & extremum_amplitude
        if (.not. trace_initialized) then
            call get_environment_variable("GRASP_TRACE_COUNT", trace_setting, status=trace_status)
            trace_enabled = trace_status == 0 .and. (trace_setting(1:1) == "1" .or. trace_setting(1:1) == "y")
            trace_initialized = .true.
        end if
        extremum_count = 1
        segment_max_amplitude = 0.0_real64
        extremum_grid_indices(1) = 1
        peak_amplitude = 0.0_real64
        do sample_index = 2, last_radial_index
            sample_amplitude = abs(radial_values(sample_index))
            absolute_sign_sum = abs(sign(1.0_real64,radial_values(sample_index))+sign(1.0_real64,radial_values(sample_index-1)))
            if (absolute_sign_sum <= working_tolerance) then
                extremum_count = extremum_count + 1
                extremum_grid_indices(extremum_count) = 0
                segment_max_amplitude = 0.0_real64
            end if
            if (sample_amplitude > segment_max_amplitude) then
                segment_max_amplitude = sample_amplitude
                extremum_grid_indices(extremum_count) = sample_index
            end if
            if (sample_amplitude <= peak_amplitude) cycle
            peak_amplitude = sample_amplitude
        end do
        extremum_index = 0
        ! Ignore segment extrema below the configured fraction of the sampled peak amplitude.
        extremum_cutoff = relative_extremum_cutoff*peak_amplitude
        do
            extremum_index = extremum_index + 1
            if (extremum_index > extremum_count) exit
            if (extremum_grid_indices(extremum_index) == 0) then
                extremum_amplitude = 0.0_real64
            else
                extremum_amplitude = abs(radial_values(extremum_grid_indices(extremum_index)))
            end if
            if (extremum_amplitude < extremum_cutoff) then
                extremum_count = extremum_count - 1
                remaining_extrema = extremum_count - extremum_index
                extremum_grid_indices(extremum_index:remaining_extrema+extremum_index) = &
                    & extremum_grid_indices(extremum_index+1:remaining_extrema+1+extremum_index)
                extremum_index = extremum_index - 1
            end if
        end do
        ! Count sign changes between retained extrema so sub-threshold oscillations do not add nodes.
        node_count = 0
        do sample_index = 2, extremum_count
            absolute_sign_sum = abs(sign(1.0_real64,radial_values(extremum_grid_indices(sample_index)))+sign(1.0_real64, &
                & radial_values(extremum_grid_indices(sample_index-1))))
            if (absolute_sign_sum > working_tolerance) cycle
            node_count = node_count + 1
        end do
        first_extremum_sign = sign(1.0_real64,radial_values(extremum_grid_indices(1)))
        if (trace_enabled .and. last_radial_index > 300 .and. count_trace_rank == 0) then
            call ensure_calculation_log()
            write (calculation_log_unit,"(A,I0,A,I0,A,ES16.8,A,ES16.8,A,I0)") "COUNT_TRACE orbital=", &
                count_context, " mtp=", last_radial_index, " emx=", peak_amplitude, " thre=", extremum_cutoff, " nodes=", &
                & node_count
            call ensure_calculation_log()
            write (calculation_log_unit,"(A,I0)") "COUNT_EXTREMA orbital=", count_context
            do sample_index = 1, extremum_count
                if (extremum_grid_indices(sample_index) > 0) then
                    call ensure_calculation_log()
                    write (calculation_log_unit,"(I0,A,I0,A,ES16.8,A,ES16.8)") &
                        sample_index, ",", extremum_grid_indices(sample_index), ",", &
                        & radial_grid(extremum_grid_indices(sample_index)), ",", &
                        & abs(radial_values(extremum_grid_indices(sample_index)))
                end if
            end do
        end if
        return
    end subroutine count_radial_nodes
end module radial_node_count
