!***********************************************************************
!                                                                      *
subroutine calculate_deferred_corrections(orbital_index)
    ! Module dependencies
    !                                                                      *
    !   Compute the deferred corrections for orbital J .                   *
    !                                                                      *
    !                                          Last updated: 18 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use atomic_calculation_state, only: dp, dq
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_index

    ! Local constants
    real(real64), parameter :: w3 = 1.0_real64/120.0_real64
    real(real64), parameter :: w2 = -15.0_real64*w3
    real(real64), parameter :: w1 = 40.0_real64*w3

    ! Local declarations
    integer :: i, mfjm3, mfjm2
    logical :: first

    ! Initial values
    !-----------------------------------------------
    !
    data first/ .true./
    !
    !   The deferred corrections for the first two points are
    !   unnecessary, because the integration always commences
    !   from the fourth point of the grid; this requires the
    !   deferred correction at the third and subsequent points
    !   only
    !
    if (first) then
        dp(:2) = 0.0_real64
        dq(:2) = 0.0_real64
        first = .false.
    end if
    !
    !   Intermediate points
    !
    mfjm3 = orbital_last_grid_index(orbital_index) - 3
    do i = 3, mfjm3
        !
        dp(i) = w3*(radial_large_components(i + 3,orbital_index) - radial_large_components(i - 2,orbital_index)) + &
            w2*(radial_large_components(i + 2,orbital_index) - radial_large_components(i - 1,orbital_index)) + &
            w1*(radial_large_components(i + 1,orbital_index) - radial_large_components(i,orbital_index))
        !
        dq(i) = w3*(radial_small_components(i + 3,orbital_index) - radial_small_components(i - 2,orbital_index)) + &
            w2*(radial_small_components(i + 2,orbital_index) - radial_small_components(i - 1,orbital_index)) + &
            w1*(radial_small_components(i + 1,orbital_index) - radial_small_components(i,orbital_index))
        !
    end do
    !
    !   Set remaining deferred corrections to zero: slopes are
    !   small in this region
    !
    mfjm2 = orbital_last_grid_index(orbital_index) - 2
    dp(mfjm2:radial_point_count) = 0.0_real64
    dq(mfjm2:radial_point_count) = 0.0_real64
    !
    return
end subroutine calculate_deferred_corrections
