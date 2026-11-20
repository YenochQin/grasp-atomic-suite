!***********************************************************************
!                                                                      *
subroutine build_variation_source_term(orbital_index)
    ! Module dependencies
    !                                                                      *
    !   This  subprogram  sets up the inhomogeneous terms for the varia-   *
    !   tion equations.                                                    *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use atomic_calculation_state, only: max_radial_points, speed_of_light_au
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    implicit none (type, external)

    ! Arguments
    integer  :: orbital_index

    ! Local declarations
    integer :: radial_index
    real(real64) :: hhc
    !-----------------------------------------------
    !
    hhc = 0.5_real64*transformed_grid_step/speed_of_light_au
    !
    !   Set up arrays TF and TG
    !
    do radial_index = 1, radial_point_count
        xu(radial_index) = -radial_small_components(radial_index,orbital_index)*hhc*radial_grid_derivative(radial_index)
        xv(radial_index) = radial_large_components(radial_index,orbital_index)*hhc*radial_grid_derivative(radial_index)
    end do
    !
    xu(:radial_point_count-1) = xu(2:radial_point_count) + xu(:radial_point_count-1)
    xv(:radial_point_count-1) = xv(2:radial_point_count) + xv(:radial_point_count-1)
    !
    return
    !
end subroutine build_variation_source_term
