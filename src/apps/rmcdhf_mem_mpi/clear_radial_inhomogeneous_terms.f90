!***********************************************************************
!                                                                      *
subroutine clear_radial_inhomogeneous_terms(orbital_index)
    ! Module dependencies
    !                                                                      *
    !   This subprogram sets the inhomogeneous terms to zero.              *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
    implicit none (type, external)

    ! Arguments
    integer  :: orbital_index
    !-----------------------------------------------
    !
    xu(:radial_point_count) = 0.0_real64
    xv(:radial_point_count) = 0.0_real64
    !
    return
end subroutine clear_radial_inhomogeneous_terms
