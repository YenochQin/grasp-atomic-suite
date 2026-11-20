!***********************************************************************
!                                                                      *
subroutine build_qed_self_energy(jstate, slfint, ucf)
    ! Module dependencies
    !                                                                      *
    !   This  routine estimates corrections to the  energy levels due to   *
    !   self-energy.                                                       *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SCREEN.      *
    !               [RCI92]: FZALF, HOVLAP.                                *
    !                                                                      *
    !                                           Last update: 30 Oct 1992   *
    !   Modified by Xinghong He                 Last update: 24 Jun 1997   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals, max_radial_points
    use atomic_calculation_state, only: max_radial_points, nuclear_charge, pi, speed_of_light_au
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
    use orbital_csf_basis, only: orbital_count, csf_count, orbital_principal_n, orbital_kappa
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use csf_quantum_queries, only: iq

    ! Procedure interfaces
    use radial_density_ratio_interface, only: radial_density_ratio
    use hydrogenic_self_energy_factor_interface, only: hydrogenic_self_energy_factor
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: jstate
    real(real64), intent(out) :: slfint(max_orbitals)
    real(real64), intent(out) :: ucf(:)

    ! Local declarations
    integer :: maxiter, orbital_index, i, csf_index, npj, kappa, last_radial_index
    real(real64) :: ptemp(max_radial_points), qtemp(max_radial_points)
    real(real64) :: ucfj, zeff, ratio, valu
    character :: npchar, nakchar*2
    !-----------------------------------------------
    !
    ! Pre-set tolerable number for iteration in finding effective
    ! nuclear charge.
    !
    maxiter = 20
    !
    ! Modified so that UCFJ describes the current eigenstate
    !
    do orbital_index = 1, orbital_count
        ucfj = 0.0_real64
        !         DO 3 I = 1,NVEC
        i = jstate
        do csf_index = 1, csf_count
            ucfj = ucfj + dble(iq(orbital_index,csf_index))*csf_mixing_coefficients(csf_index + (i - 1)*csf_count)**2
        end do
        !    3    CONTINUE
        !     print *, ucfj,'ucf'
        ucf(orbital_index) = ucfj
        ! zou    UCF(J) = UCFJ/DBLE (NCF)
    end do
    !
    do orbital_index = 1, orbital_count
        !
        npj = orbital_principal_n(orbital_index)
        !
        if (npj <= 8) then
            !
            !   Only orbitals with principal quantum number 8 or less can
            !   be treated by this section of code
            !
            kappa = orbital_kappa(orbital_index)
            !
            !   Begin by transferring the function to a temporary array
            !
            last_radial_index = orbital_last_grid_index(orbital_index)
            !
            ptemp(1) = 0.0_real64
            qtemp(1) = 0.0_real64
            do i = 2, last_radial_index
                ptemp(i) = radial_large_components(i,orbital_index)
                qtemp(i) = radial_small_components(i,orbital_index)
            end do
            !
            zeff = nuclear_charge
            ratio = radial_density_ratio(ptemp,qtemp,last_radial_index,npj,kappa,zeff)
            valu = ratio*hydrogenic_self_energy_factor(npj,kappa,zeff)/dble(npj**3)
            !
            slfint(orbital_index) = valu*zeff**4/(pi*speed_of_light_au**3)
        !        print *, 'No. orb.=',j,' Zeff = ',zeff
        !    &   , 'Scale= ',ratio
        !    &   , 'S.E. = ',slfint(j)*2*13.6058,slfint(j)/ratio*2*13.6058
        !
        else
            !
            !   The self-energy for orbitals with principal quantum number
            !   greater than 8 is set to zero
            !
            slfint(orbital_index) = 0.0_real64
            !
        end if
        !
    end do
    !
    !   Deallocate storage for the `generalised occupation numbers'
    !
    !
    return
end subroutine build_qed_self_energy


!***********************************************************************
!                                                                      *
real(real64) function radial_density_ratio (large_component_values, small_component_values, last_radial_index, &
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
    !
    !...Switches:
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
    integer , intent(in) :: last_radial_index
    integer  :: principal_n
    integer  :: kappa
    real(real64)  :: effective_nuclear_charge
    real(real64), intent(in) :: large_component_values(max_radial_points)
    real(real64), intent(in) :: small_component_values(max_radial_points)

    ! Local declarations
    integer :: mtph, i, k
    real(real64) :: eh, pzh, result, result1
    !-----------------------------------------------
    !
    !
    !
    !   Set up the hydrogenic orbital
    !
    call dcbsrw (principal_n, kappa, effective_nuclear_charge, eh, pzh, ph, qh, mtph)
    !
    !   Compute the overlap
    !
    mtp = min(mtph,last_radial_index)
    do i = 2, mtp
        if (radial_grid_derivative(i) > 0.0219_real64) cycle
        k = i
    end do
    mtp = k
    ta(1) = 0.0_real64
    ta(2:mtp) = (large_component_values(2:mtp)*large_component_values(2:mtp)+small_component_values(2:mtp)* &
        & small_component_values(2:mtp))*radial_grid_derivative(2:mtp)
    !         TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
    call quad (result)
    ta(2:mtp) = (ph(2:mtp)*ph(2:mtp)+qh(2:mtp)*qh(2:mtp))*radial_grid_derivative(2:mtp)
    call quad (result1)
    !
    radial_density_ratio = result/result1
    !
    return
end function radial_density_ratio
