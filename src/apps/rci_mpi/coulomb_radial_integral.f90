!***********************************************************************
!                                                                      *
real(real64) function coulomb_radial_integral (recoupling_coefficient, orbital_a, orbital_c, rbd, orbital_b, orbital_d, &
    & multipole_rank, iw)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This routine evaluates the transverse interaction integrals.  If   *
    !   IW = 0, it calulates the U(r1,r2) integral; otherwise, it calcu-   *
    !   lates  R bar (k; a c | b d ; w)  with  w = wac  if  IW = 1,  and   *
    !   w = wbd  if IW = 2.                                                *
    !                                                                      *
    !   Call(s) to: [LIB92]: QUAD.                                         *
    !               [RCI92]: ZKF.                                          *
    !                                           Last update: 15 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_radial_points
    use bessel_function_workspace, only: bessj, bessn, bj, max_radial_points, tc, td, wij
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use radial_grid_state, only: radial_derivative_over_radius
    use orbital_csf_basis, only: max_orbitals, orbital_principal_n, orbital_symmetry_labels
    use radial_quadrature_workspace,  only: mtp, ta, tb
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii

    ! Procedure interfaces
    use evaluate_breit_hartree_z_interface, only: evaluate_breit_hartree_z

    ! Module dependencies
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_a, orbital_c
    integer, intent(in) :: orbital_b, orbital_d
    integer, intent(in) :: multipole_rank
    integer, intent(in) :: iw
    real(real64), intent(in) :: recoupling_coefficient(max_radial_points)
    real(real64), intent(in) :: rbd(max_radial_points)

    ! Local declarations
    integer :: mxrbd, mxrac, i
    real(real64) :: result
    !-----------------------------------------------
    !
    mxrbd = min (orbital_last_grid_index(orbital_b),orbital_last_grid_index(orbital_d))
    mxrac = min (orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
    !
    if (iw == 0) then
        !
        !   IW = 0
        !
        ta(:mxrac) = recoupling_coefficient(:mxrac)
        mtp = mxrac
        call evaluate_breit_hartree_z (multipole_rank, orbital_a, orbital_c)
        mtp = min(mtp,mxrbd)
        ta(1) = 0.0_real64
        ta(2:mtp) = rbd(2:mtp)*tb(2:mtp)*radial_derivative_over_radius(2:mtp)
    !
    else
        !
        !   IW = 1,2
        !
        ta(:mxrac) = recoupling_coefficient(:mxrac)*(1.0_real64 + bessj(1,iw,:mxrac))
        mtp = mxrac
        call evaluate_breit_hartree_z (multipole_rank, orbital_a, orbital_c)
        mtp = min(mtp,mxrbd)
        ta(1) = 0.0_real64
        ta(2:mtp) = rbd(2:mtp)*(1.0_real64 + bessn(1,iw,2:mtp))*tb(2:mtp)*radial_derivative_over_radius(2:&
            mtp)
        !
    end if
    !
    call quad (result)
    coulomb_radial_integral = result
    !
    !   Debug printout if option set
    !
    if (ldbpr(11)) then
        call ensure_debug_output()
        write (angular_debug_unit, '("_ (",1i2,")"/,"R     (",1i2,1a2,",",1i2,1a2,"|",1i2,1a2,",",1i2,1a2,";",1i2,&
            &") = ",1p,d19.12)') multipole_rank, orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), &
            & orbital_principal_n(orbital_c), &
            & orbital_symmetry_labels(orbital_c), orbital_principal_n(orbital_b)&
            , orbital_symmetry_labels(orbital_b), orbital_principal_n(orbital_d), orbital_symmetry_labels(orbital_d), iw, result
    end if
    !
    return
    !
    return
    !
end function coulomb_radial_integral
