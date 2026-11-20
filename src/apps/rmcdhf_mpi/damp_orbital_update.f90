!***********************************************************************
!                                                                      *
subroutine damp_orbital_update(orbital_index, orbital_sign_reversal_count, orbital_damping_factor)
    ! Module dependencies
    !                                                                      *
    !   This subroutine damps the orbital wave function with index J. it   *
    !   also  stores  the  previous  determination  of this orbital.       *
    !                                                                      *
    !   Call(s) to: [LIB92]: COUNT, RINT.                                  *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_radial_points
    use scf_damping_factors, only: max_orbitals, mixing_damping_factors, orbital_damping_factors
    use atomic_calculation_state, only: accy, max_radial_points
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
    use orbital_inversion_flags, only: max_orbitals, noinvt
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use weighted_radial_overlap, only: rint
    !      USE count_I
    use radial_node_count, only: count
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_index
    integer, intent(inout) :: orbital_sign_reversal_count
    real(real64), intent(in) :: orbital_damping_factor

    ! Local declarations
    integer :: mtpo, mtpn, mtp, i, nncff, last_radial_index
    real(real64) :: eps, factor, poldi, qoldi, original_norm_squared, normalization_factor, sgn
    logical :: check
    !-----------------------------------------------
    !
    !
    !   Initialization
    !
    !ww      EPS = 0.1D 00*ACCY
    eps = 0.01_real64*accy
    check = .not.noinvt(orbital_index)
    !
    !   Damp orbital J using the damping factor  ABS (ODAMP(J)):  ODAMP(J)
    !   is restricted to the open interval (-1,1) by  DATSCF ; the meaning
    !   of a  negative  ODAMP(J)  is that  ABS (ODAMP(J))  is the constant
    !   damping factor;  a positive  ODAMP(J)  can be reset by  subroutine
    !   dampck
    !
    !   Store previous determination of this orbital in arrays P, Q
    !
    !XHH odampj goes to the argument
    !      ODAMPJ = ABS (ODAMP(J))
    !
    if (orbital_damping_factor > eps) then
        !
        factor = 1.0_real64 - orbital_damping_factor
        !
        orbital_origin_parameters(orbital_index) = factor*p0 + orbital_damping_factor*orbital_origin_parameters(orbital_index)
        !
        mtpo = orbital_last_grid_index(orbital_index)
        mtpn = mtp0
        mtp0 = mtpo
        !
        mtp = max(mtpn,mtpo)
        do i = 1, mtp
            poldi = radial_large_components(i,orbital_index)
            radial_large_components(i,orbital_index) = factor*p(i) + orbital_damping_factor*radial_large_components(i,orbital_index)
            p(i) = poldi
            qoldi = radial_small_components(i,orbital_index)
            radial_small_components(i,orbital_index) = factor*q(i) + orbital_damping_factor*radial_small_components(i,orbital_index)
            q(i) = qoldi
        end do
        !
        !   Compute normalization factor
        !
        orbital_last_grid_index(orbital_index) = mtp
        original_norm_squared = rint(orbital_index,orbital_index,0)
        normalization_factor = 1.0_real64/sqrt(original_norm_squared)
        !
        !   Determine if inversion is necessary
        !
        if (check) then
            call count (radial_large_components(:max_radial_points,orbital_index), mtp, nncff, sgn)
            if (sgn < 0.0_real64) then
                orbital_sign_reversal_count = orbital_sign_reversal_count + 1
                normalization_factor = -normalization_factor
            end if
        end if
        !
        orbital_origin_parameters(orbital_index) = orbital_origin_parameters(orbital_index)*normalization_factor
        radial_large_components(:mtp,orbital_index) = normalization_factor*radial_large_components(:mtp,orbital_index)
        radial_small_components(:mtp,orbital_index) = normalization_factor*radial_small_components(:mtp,orbital_index)
        !
        !   Find new MF(J)
        !
        last_radial_index = mtp + 1
        do
            last_radial_index = last_radial_index - 1
            if (abs(radial_large_components(last_radial_index,orbital_index)) < eps) then
                radial_large_components(last_radial_index,orbital_index) = 0.0_real64
                radial_small_components(last_radial_index,orbital_index) = 0.0_real64
            else
                orbital_last_grid_index(orbital_index) = last_radial_index
                exit
            end if
        end do
    !
    else
        !
        orbital_origin_parameters(orbital_index) = p0
        !
        mtpo = orbital_last_grid_index(orbital_index)
        mtpn = mtp0
        mtp0 = mtpo
        !
        mtp = max(mtpn,mtpo)
        do i = 1, mtp
            poldi = radial_large_components(i,orbital_index)
            radial_large_components(i,orbital_index) = p(i)
            p(i) = poldi
            qoldi = radial_small_components(i,orbital_index)
            radial_small_components(i,orbital_index) = q(i)
            q(i) = qoldi
        end do
        !
        orbital_last_grid_index(orbital_index) = mtp
        !
    end if
    !
    return
end subroutine damp_orbital_update
