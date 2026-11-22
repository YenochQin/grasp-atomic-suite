! Evaluate Hamiltonian corrections that are deferred until after diagonalization.
module calculate_deferred_corrections_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_deferred_corrections
contains

    !   Compute the deferred corrections for orbital J .
    !
    !                                          Last updated: 18 Dec 1992
    subroutine calculate_deferred_corrections(orbital_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: large_deferred_correction, small_deferred_correction
        use radial_grid_state, only: radial_point_count
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index

        ! Local constants
        real(real64), parameter :: outer_stencil_weight = 1.0_real64/120.0_real64
        real(real64), parameter :: inner_stencil_weight = -15.0_real64*outer_stencil_weight
        real(real64), parameter :: central_stencil_weight = 40.0_real64*outer_stencil_weight

        ! Local declarations
        integer :: radial_index, last_orbital_grid_index_minus_three, last_orbital_grid_index_minus_two
        logical :: first

        ! Initial values
        data first/ .true./
        ! The deferred corrections for the first two points are
        ! unnecessary, because the integration always commences
        ! from the fourth point of the grid; this requires the
        ! deferred correction at the third and subsequent points
        ! only
        if (first) then
            large_deferred_correction(:2) = 0.0_real64
            small_deferred_correction(:2) = 0.0_real64
            first = .false.
        end if
        ! Intermediate points
        last_orbital_grid_index_minus_three = orbital_last_grid_index(orbital_index) - 3
        do radial_index = 3, last_orbital_grid_index_minus_three
            large_deferred_correction(radial_index) = outer_stencil_weight*(radial_large_components(radial_index + 3, &
                orbital_index) - &
                & radial_large_components(radial_index - 2,orbital_index)) + &
                inner_stencil_weight*(radial_large_components(radial_index + 2,orbital_index) - &
                    & radial_large_components(radial_index - 1,orbital_index)) + &
                central_stencil_weight*(radial_large_components(radial_index + 1,orbital_index) - &
                    & radial_large_components(radial_index,orbital_index))
            small_deferred_correction(radial_index) = outer_stencil_weight*(radial_small_components(radial_index + 3, &
                orbital_index) - &
                & radial_small_components(radial_index - 2,orbital_index)) + &
                inner_stencil_weight*(radial_small_components(radial_index + 2,orbital_index) - &
                    & radial_small_components(radial_index - 1,orbital_index)) + &
                central_stencil_weight*(radial_small_components(radial_index + 1,orbital_index) - &
                    & radial_small_components(radial_index,orbital_index))
        end do
        ! Set remaining deferred corrections to zero: slopes are
        ! small in this region
        last_orbital_grid_index_minus_two = orbital_last_grid_index(orbital_index) - 2
        large_deferred_correction(last_orbital_grid_index_minus_two:radial_point_count) = 0.0_real64
        small_deferred_correction(last_orbital_grid_index_minus_two:radial_point_count) = 0.0_real64
        return
    end subroutine calculate_deferred_corrections
end module calculate_deferred_corrections_interface
