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
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_index

    ! Local declarations
    integer :: mtp, i
    real(real64) :: scmea, deltao
    !-----------------------------------------------
    !
    !
    scmea = 0.0_real64
    mtp = min(mtp0,orbital_last_grid_index(orbital_index))
    do i = 1, mtp
        deltao = abs(p(i)-radial_large_components(i,orbital_index)) + abs(q(i)-radial_small_components(i,orbital_index))
        scmea = dmax1(deltao,scmea)
    end do
    scnsty(orbital_index) = scmea*sqrt(ucf(orbital_index))
    !
    return
end subroutine evaluate_radial_self_consistency
