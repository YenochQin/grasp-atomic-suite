module damp_orbital_update_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: damp_orbital_update
contains

    !***********************************************************************
    !                                                                      *
    subroutine damp_orbital_update(orbital_index, orbital_sign_reversal_count, orbital_damping_factor)
        ! Module dependencies
        !                                                                      *
        !   This subroutine damps the orbital wave function with index J. it   *
        !   also  stores  the  previous  determination  of this orbital.       *
        !                                                                      *
        !   Call(s) to: [LIB92]: COUNT, RINT.                                  *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_radial_points
        use atomic_calculation_state, only: working_tolerance, max_radial_points
        use radial_dirac_workspace, only: max_radial_points, trial_orbital_last_grid_index, trial_large_component, &
            trial_large_origin_parameter, trial_small_component
        use orbital_inversion_flags, only: suppress_radial_inversion
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap
        !      USE count_I
        use radial_node_count, only: count_radial_nodes
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        integer, intent(inout) :: orbital_sign_reversal_count
        real(real64), intent(in) :: orbital_damping_factor

        ! Local declarations
        integer :: old_last_grid_index, trial_last_grid_index, common_last_grid_index, radial_index, radial_node_count, &
            & last_radial_index
        real(real64) :: damping_threshold, trial_weight, old_large_component, old_small_component, original_norm_squared, &
            & normalization_factor, first_extremum_sign
        logical :: allow_sign_reversal
        !-----------------------------------------------
        !
        !
        !   Initialization
        !
        !ww      EPS = 0.1D 00*ACCY
        damping_threshold = 0.01_real64*working_tolerance
        allow_sign_reversal = .not.suppress_radial_inversion(orbital_index)
        !
        !   Damp orbital J using the damping factor  ABS (ODAMP(J)):  ODAMP(J)
        !   is restricted to the open interval (-1,1) by  DATSCF ; the meaning
        !   of a  negative  ODAMP(J)  is that  ABS (ODAMP(J))  is the constant
        !   damping factor;  a positive  ODAMP(J)  can be reset by  subroutine
        !   dampck
        !
        !   Store previous determination of this orbital in arrays P, Q
        !
        !XHH odampj goes to the argument
        !      ODAMPJ = ABS (ODAMP(J))
        !
        if (orbital_damping_factor > damping_threshold) then
            !
            trial_weight = 1.0_real64 - orbital_damping_factor
            !
            orbital_origin_parameters(orbital_index) = trial_weight*trial_large_origin_parameter + orbital_damping_factor* &
                orbital_origin_parameters(orbital_index)
            !
            old_last_grid_index = orbital_last_grid_index(orbital_index)
            trial_last_grid_index = trial_orbital_last_grid_index
            trial_orbital_last_grid_index = old_last_grid_index
            !
            common_last_grid_index = max(trial_last_grid_index,old_last_grid_index)
            do radial_index = 1, common_last_grid_index
                old_large_component = radial_large_components(radial_index,orbital_index)
                radial_large_components(radial_index,orbital_index) = trial_weight*trial_large_component(radial_index) + &
                    & orbital_damping_factor* &
                    radial_large_components(radial_index, &
                    & orbital_index)
                trial_large_component(radial_index) = old_large_component
                old_small_component = radial_small_components(radial_index,orbital_index)
                radial_small_components(radial_index,orbital_index) = trial_weight*trial_small_component(radial_index) + &
                    & orbital_damping_factor* &
                    radial_small_components(radial_index, &
                    & orbital_index)
                trial_small_component(radial_index) = old_small_component
            end do
            !
            !   Compute normalization factor
            !
            orbital_last_grid_index(orbital_index) = common_last_grid_index
            original_norm_squared = integrate_weighted_radial_overlap(orbital_index,orbital_index,0)
            normalization_factor = 1.0_real64/sqrt(original_norm_squared)
            !
            !   Determine if inversion is necessary
            !
            if (allow_sign_reversal) then
                call count_radial_nodes (radial_large_components(:max_radial_points,orbital_index), common_last_grid_index, &
                    & radial_node_count, first_extremum_sign)
                if (first_extremum_sign < 0.0_real64) then
                    orbital_sign_reversal_count = orbital_sign_reversal_count + 1
                    normalization_factor = -normalization_factor
                end if
            end if
            !
            orbital_origin_parameters(orbital_index) = orbital_origin_parameters(orbital_index)*normalization_factor
            radial_large_components(:common_last_grid_index,orbital_index) = normalization_factor* &
                & radial_large_components(:common_last_grid_index,orbital_index)
            radial_small_components(:common_last_grid_index,orbital_index) = normalization_factor* &
                & radial_small_components(:common_last_grid_index,orbital_index)
            !
            !   Find new MF(J)
            !
            last_radial_index = common_last_grid_index + 1
            do
                last_radial_index = last_radial_index - 1
                if (abs(radial_large_components(last_radial_index,orbital_index)) < damping_threshold) then
                    radial_large_components(last_radial_index,orbital_index) = 0.0_real64
                    radial_small_components(last_radial_index,orbital_index) = 0.0_real64
                else
                    orbital_last_grid_index(orbital_index) = last_radial_index
                    exit
                end if
            end do
        !
        else
            !
            orbital_origin_parameters(orbital_index) = trial_large_origin_parameter
            !
            old_last_grid_index = orbital_last_grid_index(orbital_index)
            trial_last_grid_index = trial_orbital_last_grid_index
            trial_orbital_last_grid_index = old_last_grid_index
            !
            common_last_grid_index = max(trial_last_grid_index,old_last_grid_index)
            do radial_index = 1, common_last_grid_index
                old_large_component = radial_large_components(radial_index,orbital_index)
                radial_large_components(radial_index,orbital_index) = trial_large_component(radial_index)
                trial_large_component(radial_index) = old_large_component
                old_small_component = radial_small_components(radial_index,orbital_index)
                radial_small_components(radial_index,orbital_index) = trial_small_component(radial_index)
                trial_small_component(radial_index) = old_small_component
            end do
            !
            orbital_last_grid_index(orbital_index) = common_last_grid_index
            !
        end if
        !
        return
    end subroutine damp_orbital_update
end module damp_orbital_update_interface
