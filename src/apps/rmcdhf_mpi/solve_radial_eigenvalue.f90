!***********************************************************************
!                                                                      *
real(real64) function solve_radial_eigenvalue (orbital_index)
    ! Module dependencies
    !                                                                      *
    !   This function computes an estimate of the energy of orbital J .    *
    !                                                                      *
    !   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use atomic_calculation_state, only: max_radial_points, speed_of_light_au
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use orbital_csf_basis, only: max_orbitals, orbital_kappa
    use scf_radial_potentials, only: max_radial_points, xp, xq, yp
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use radial_quadrature, only: quad
    use radial_wavefunction_derivatives, only: dpbdt
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_index

    ! Local declarations
    integer :: i
    real(real64) :: piece1, piece2, piece3, piece4, piece5
    !-----------------------------------------------
    !
    !
    !   Initialization
    !
    mtp = orbital_last_grid_index(orbital_index)
    !
    !   Exchange term
    !
    ta(1) = 0.0_real64
    do i = 2, mtp
        ta(i) = (radial_large_components(i,orbital_index)*xq(i)-radial_small_components(i, &
            & orbital_index)*xp(i))*radial_derivative_over_radius(i)
    end do
    call quad (piece1)
    piece1 = speed_of_light_au*piece1
    !
    !   Direct term
    !
    ta(1) = 0.0_real64
    do i = 2, mtp
        ta(i) = (radial_large_components(i,orbital_index)**2 + radial_small_components(i, &
            & orbital_index)**2)*yp(i)*radial_derivative_over_radius(i)
    end do
    call quad (piece2)
    !
    !   Kinetic energy terms
    !
    ta(1) = 0.0_real64
    do i = 2, mtp
        ta(i) = radial_small_components(i,orbital_index)**2*radial_grid_derivative(i)
    end do
    call quad (piece3)
    piece3 = 2.0_real64*speed_of_light_au*speed_of_light_au*piece3
    !
    ta(1) = 0.0_real64
    do i = 2, mtp
        ta(i) = (radial_large_components(i,orbital_index)*radial_small_components(i, &
            & orbital_index))*radial_derivative_over_radius(i)
    end do
    call quad (piece4)
    piece4 = -2.0_real64*dble(orbital_kappa(orbital_index))*speed_of_light_au*piece4
    !
    call dpbdt (orbital_index)
    ta(1) = 0.0_real64
    do i = 2, mtp
        ta(i) = radial_large_components(i,orbital_index)*tb(i) - radial_small_components(i,orbital_index)*ta(i)
    end do
    call quad (piece5)
    piece5 = speed_of_light_au*piece5/transformed_grid_step
    !
    !   Assembly
    !
    solve_radial_eigenvalue = piece1 + piece2 + piece3 + piece4 + piece5
    !
    return
    !
end function solve_radial_eigenvalue
