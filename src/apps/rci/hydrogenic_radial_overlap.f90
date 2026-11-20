!***********************************************************************
!                                                                      *
real(real64) function hydrogenic_radial_overlap (large_component_values, small_component_values, last_radial_index, &
    & principal_n, kappa, &
    & effective_nuclear_charge)
    ! Module dependencies
    !                                                                      *
    !   This subprogram computes the overlap of the orbital tabulated in   *
    !   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
    !   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
    !                                                                      *
    !   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_radial_points
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use hydrogenic_radial_components, only: ph, qh
    use radial_quadrature_workspace, only: mtp, ta
    use dirac_coulomb_radial_orbital, only: dcbsrw
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: last_radial_index
    integer  :: principal_n
    integer  :: kappa
    real(real64)  :: effective_nuclear_charge
    real(real64), intent(in) :: large_component_values(max_radial_points)
    real(real64), intent(in) :: small_component_values(max_radial_points)

    ! Local declarations
    integer :: mtph, i
    real(real64) :: eh, pzh, result
    !-----------------------------------------------
    !
    !   Set up the hydrogenic orbital
    !
    call dcbsrw (principal_n, kappa, effective_nuclear_charge, eh, pzh, ph, qh, mtph)
    !
    !   Compute the overlap
    !
    mtp = min(mtph,last_radial_index)
    ta(1) = 0.0_real64
    ta(2:mtp) = (large_component_values(2:mtp)*ph(2:mtp)+small_component_values(2:mtp)*qh(2:mtp))*radial_grid_derivative(2:mtp)
    call quad (result)
    !
    hydrogenic_radial_overlap = result
    !
    return
end function hydrogenic_radial_overlap
