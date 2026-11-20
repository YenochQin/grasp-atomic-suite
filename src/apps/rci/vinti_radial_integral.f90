!***********************************************************************
!                                                                      *
real(real64) function vinti_radial_integral (first_orbital, second_orbital)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   The value of this  function is the one-electron integral V (J,K)   *
    !   for  orbitals  J, K. The analytical expression for this quantity   *
    !   is given as  eq (3.23) in  F A  Parpia, M. Tong and C F Fischer,   *
    !   to appear.                                                         *
    !                                                                      *
    !   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
    !                                                                      *
    !   Written by M Tong and F A Parpia,     Last revision: 15 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use orbital_csf_basis, only: max_orbitals, orbital_kappa, orbital_principal_n, orbital_symmetry_labels
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use radial_wavefunction_derivatives, only: dpbdt
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: first_orbital
    integer, intent(in) :: second_orbital

    ! Local declarations
    integer :: i, kpj, kpk, ifact1, ifact2
    real(real64) :: piece1, fact1, fact2, piece2
    !-----------------------------------------------
    !
    mtp = max(orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
    !
    !   Piece involving derivatives
    !
    call dpbdt (second_orbital)
    ta(1) = 0.0_real64
    do i = 2, mtp
        ta(i) = radial_large_components(i,first_orbital)*ta(i) + radial_small_components(i,first_orbital)*tb(i)
    end do
    call quad (piece1)
    piece1 = piece1/transformed_grid_step
    !
    !   Pieces not involving derivatives
    !
    kpj = orbital_kappa(first_orbital)
    kpk = orbital_kappa(second_orbital)
    ifact1 = kpj*(kpj + 1) - kpk*(kpk + 1)
    fact1 = 0.5_real64*dble(ifact1)
    ifact2 = (-kpj*((-kpj) + 1)) + kpk*((-kpk) + 1)
    fact2 = 0.5_real64*dble(ifact2)
    ta(1) = 0.0_real64
    do i = 2, mtp
        ta(i) = radial_derivative_over_radius(i)*(fact1*radial_large_components(i,first_orbital)*radial_large_components(i, &
            & second_orbital) + &
            & fact2*radial_small_components(i,first_orbital)*radial_small_components(i,second_orbital))
    end do
    call quad (piece2)
    !
    vinti_radial_integral = piece1 - piece2
    !
    !   Debug printout
    !
    if (ldbpr(6)) then
        call ensure_debug_output()
        write (angular_debug_unit, '(/,"VINTI: V (",1i2,1a2,",",1i2,1a2,") = ",1p,d19.12)') &
            orbital_principal_n(first_orbital), orbital_symmetry_labels(first_orbital), orbital_principal_n(second_orbital), &
            & orbital_symmetry_labels(second_orbital), vinti_radial_integral
    end if
    !
    return
    !
    return
    !
end function vinti_radial_integral
