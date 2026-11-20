module evaluate_radial_self_consistency_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_radial_self_consistency
contains

    !***********************************************************************
    !                                                                      *
    subroutine evaluate_radial_self_consistency(orbital_index)
        ! Module dependencies
        !                                                                      *
        !   This routine computes the weighted self-consistency of orbital J   *
        !                                                                      *
        !   Written by Farid A Parpia, at OXFORD    Last update: 08 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_dirac_workspace, only: trial_orbital_last_grid_index, trial_large_component, trial_small_component
        use scf_iteration_state, only: orbital_occupation_weights, orbital_self_consistency_errors
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index

        ! Local declarations
        integer :: common_last_grid_index, radial_index
        real(real64) :: maximum_component_difference, component_difference
        !-----------------------------------------------
        !
        !
        maximum_component_difference = 0.0_real64
        common_last_grid_index = min(trial_orbital_last_grid_index,orbital_last_grid_index(orbital_index))
        do radial_index = 1, common_last_grid_index
            component_difference = abs(trial_large_component(radial_index)-radial_large_components(radial_index,orbital_index)) + &
                & abs(trial_small_component(radial_index)-radial_small_components(radial_index,orbital_index))
            maximum_component_difference = dmax1(component_difference,maximum_component_difference)
        end do
        orbital_self_consistency_errors(orbital_index) = maximum_component_difference* &
            sqrt(orbital_occupation_weights(orbital_index))
        !
        return
    end subroutine evaluate_radial_self_consistency
end module evaluate_radial_self_consistency_interface
