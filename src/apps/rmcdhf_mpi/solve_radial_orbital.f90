!***********************************************************************
!                                                                      *
subroutine solve_radial_orbital(orbital_index, solve_failed, orbital_sign_reversal_count, matching_grid_index, nnp)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This subroutine  performs  step  2 in Algorithm 5.2 and 5.3 of C   *
    !   Froese Fischer, Comput Phys Rep 3 (1986) 295. Some minor changes   *
    !   have been made.                                                    *
    !                                                                      *
    !   Arguments:                                                         *
    !                                                                      *
    !      J     : (Input) The serial number of the orbital                *
    !      JP    : (Output) The join point                                 *
    !      FAIL  : (Output) If .true., the iterations did not yield an     *
    !              acceptable solution (methods 1 and 2)                   *
    !                                                                      *
    !   Call(s) to: [RSCF92]: COUNT, ESTIM, IN, NEWE, OUT, cofpot,         *
    !                         SETXUV, SETXV, SETXZ, START                  *
    !               [LIB92]: DCBSRW, QUAD, SETPOT.                         *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_radial_points
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use atomic_calculation_state, only: speed_of_light_au, nsolv
    use radial_grid_state, only: radial_grid_derivative
    use radial_dirac_workspace, only: p, q, p0, q0, mtp0
    use orbital_inversion_flags, only: max_orbitals, noinvt
    use expected_orbital_nodes, only: max_orbitals, nnodep
    use orbital_csf_basis, only: max_orbitals, orbital_energy_parameters, orbital_kappa, orbital_l, orbital_principal_n, &
        orbital_symmetry_labels
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use scf_radial_potentials, only: max_radial_points, xp, xq, yp

    ! Procedure interfaces
    use estimate_orbital_energy_interface, only: estimate_orbital_energy
    use solve_radial_eigenvalue_interface, only: solve_radial_eigenvalue

    ! Module dependencies
    use dirac_coulomb_radial_orbital, only: dcbsrw
    use radial_dirac_coefficients, only: setpot

    ! Procedure interfaces
    use clear_radial_inhomogeneous_terms_interface, only: clear_radial_inhomogeneous_terms

    ! Module dependencies
    use radial_origin_solution, only: start

    ! Procedure interfaces
    use integrate_radial_outward_interface, only: integrate_radial_outward
    use integrate_radial_inward_interface, only: integrate_radial_inward
    use build_radial_inhomogeneous_terms_interface, only: build_radial_inhomogeneous_terms

    ! Module dependencies
    use radial_quadrature, only: quad

    ! Procedure interfaces
    use build_variation_source_term_interface, only: build_variation_source_term
    use print_radial_wavefunction_interface, only: print_radial_wavefunction

    ! Module dependencies
    use radial_node_count, only: count

    ! Procedure interfaces
    use update_trial_orbital_energy_interface, only: update_trial_orbital_energy
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_index
    integer, intent(inout) :: orbital_sign_reversal_count
    integer, intent(out) :: matching_grid_index
    integer, intent(out) :: nnp
    logical  :: solve_failed

    ! Local declarations
    integer :: kount, nakj, nloops, radial_equation_type, mtph, mtpi, i, mtpc, mtpv, mtpvc, &
        loc1, loc2, loc, node_count_difference, effective_principal_n
    real(real64) :: large_component_values(max_radial_points), small_component_values(max_radial_points), &
        & pv(max_radial_points), qv(max_radial_points)
    real(real64) :: twocsq, elast, tenej, energy_correction, first_oscillation_sign, qjpoh, qjpih, qjpoi, &
        qjpii, dnorm, alfa, p0h, origin_slope_parameter, small_component_origin_coefficient, qjpov, qjpiv, crnorm, dvnorm, aa, &
        bb, cc, discr, qq, root1, root2, pmx, apm, pati, abpi, ratio, test1, &
        test2, dele
    logical :: check
    !-----------------------------------------------
    !
    !   Initialization
    !
    check = .not.noinvt(orbital_index)
    solve_failed = .false.
    kount = 0
    nakj = orbital_kappa(orbital_index)
    !     TWOCSQ = 2.0D 00*137.036**2
    twocsq = 2.0_real64*speed_of_light_au*speed_of_light_au
    !
    !   Debug header
    !
    if (ldbpr(22)) then
        call ensure_debug_output()
        write (angular_debug_unit, '(/,/," Debug printout active; orbital: ",1i2,1a2)') &
            orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index)
    end if
    !
    nloops = max(nsolv,3*orbital_principal_n(orbital_index))
    !
    call estimate_orbital_energy (orbital_index)
    elast = orbital_energy_parameters(orbital_index)
    orbital_energy_parameters(orbital_index) = solve_radial_eigenvalue(orbital_index)
    !
    !   Checks on lower bounds
    !
    if (orbital_energy_parameters(orbital_index) < epsmin) then
        if (orbital_energy_parameters(orbital_index) > 0.0_real64) then
            epsmin = orbital_energy_parameters(orbital_index)
        else
            write (*, '(" E(",1i2,1a2,") = ",1p,d11.4,"; adjusted to EPSMIN")') &
                orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), orbital_energy_parameters(orbital_index)
            orbital_energy_parameters(orbital_index) = epsmin
            if (abs(emin - elast)<=1.0e-6_real64 .and. method(orbital_index)<=2) then
                call dcbsrw (orbital_principal_n(orbital_index), nakj, zinf, orbital_energy_parameters(orbital_index), p0, p, q, &
                    & mtp0)
                write (*, '(" Returned hydrogenic function for ",1i2,1a2," with"," effective charge ",f7.3)') &
                    orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), zinf
                return
            end if
        end if
    end if
    !
    !   Check on upper bound
    !
    if (method(orbital_index) <= 2) then
        if (orbital_energy_parameters(orbital_index) > epsmax) then
            tenej = 10.0_real64*orbital_energy_parameters(orbital_index)
            epsmax = min(tenej,twocsq)
            emax = epsmax
            if (orbital_energy_parameters(orbital_index) > twocsq) then
                write (*, '(" E(",1i2,1a2,") = ",1p,d11.4,"; adjusted to TWOCSQ")') &
                    orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                    & orbital_energy_parameters(orbital_index)
                orbital_energy_parameters(orbital_index) = twocsq
            end if
        end if
    end if
    !
    !   Iteration loop begins here
    !
    do
        kount = kount + 1
        !
        p0 = orbital_origin_parameters(orbital_index)
        !
        if (kount > 1) then
            !
            !   Check that bounds are ordered correctly
            !
            if (epsmax <= epsmin) write (*, '(" Warning: difficulty with node-counting procedure"/,&
                &" lower bound on energy (",1p,d11.4,") exceeds upper"," bound (",1d11.4,"; E(",1i2,1a2,") = ",1d11.4)') &
                epsmin, epsmax, orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                & orbital_energy_parameters(orbital_index&
                )
            !
            if (kount>nloops .or. epsmax-epsmin<1.0_real64/dble(orbital_principal_n(orbital_index))**3) then
                write (*, '(" Method ",1i1," unable to solve for ",1i2,1a2," orbital")') &
                    method(orbital_index), orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index)
                write (*, '(" Iteration number: ",1i2,", limit: ",1i2,/," Present estimate of P0; ",1d21.14,/,&
                    &" Present estimate of E(J): ",1d21.14,", DELEPS: ",1d21.14,/," Lower bound on energy: ",1d21.14,&
                    &", upper bound: ",1d21.14,/," Join point: ",1i4,", Maximum tabulation point:",1i4,/,&
                    &" Number of nodes counted: ",1i2,", Correct number: ",1i2,/," Sign of P at first oscillation: ",&
                    &f3.0)') kount - 1, nloops, p0, orbital_energy_parameters(orbital_index), energy_correction, epsmin, epsmax&
                    , matching_grid_index, mtp, nnp, nnodep(orbital_index), first_oscillation_sign
                solve_failed = .true.
                return
            end if
        end if
        !
        !   Set up arrays TF and TG; find join point
        !
        call setpot (orbital_index, matching_grid_index)
        !
        !   Set right-hand side to zero to form homogeneous equations;
        !   integrate homogeneous equations outwards and inwards; store
        !   small component at join point each time
        !
        call clear_radial_inhomogeneous_terms (orbital_index)
        radial_equation_type = 1
        call start (orbital_index, radial_equation_type, p0, large_component_values, q0, small_component_values)
        call integrate_radial_outward (orbital_index, matching_grid_index, large_component_values, small_component_values)
        qjpoh = small_component_values(matching_grid_index)
        call integrate_radial_inward (orbital_index, matching_grid_index, large_component_values, small_component_values, mtph)
        qjpih = small_component_values(matching_grid_index)
        !
        !   Set up right-hand side for inhomogeneous equations; integrate
        !   inhomogeneous equations outwards and inwards; store small
        !   component at join point each time
        !
        call build_radial_inhomogeneous_terms (orbital_index)
        radial_equation_type = 2
        call start (orbital_index, radial_equation_type, p0, p, q0, q)
        call integrate_radial_outward (orbital_index, matching_grid_index, p, q)
        qjpoi = q(matching_grid_index)
        call integrate_radial_inward (orbital_index, matching_grid_index, p, q, mtpi)
        qjpii = q(matching_grid_index)
        !
        !   Determine energy adjustment for methods 1 and 2
        !
        if (method(orbital_index) <= 2) then
            ta(1) = 0.0_real64
            ta(2:mtpi) = (p(2:mtpi)**2+q(2:mtpi)**2)*radial_grid_derivative(2:mtpi)
            mtp = mtpi
            call quad (dnorm)
            !        DELEPS = 137.036*P(JP)*(QJPII-QJPOI)/DNORM
            energy_correction = speed_of_light_au*p(matching_grid_index)*(qjpii - qjpoi)/dnorm
        end if
        !
        !   Generate the continuous solution
        !
        mtpc = max(mtph,mtpi)
        alfa = -(qjpii - qjpoi)/(qjpih - qjpoh)
        p0h = p0
        p0 = p0*(1.0_real64 + alfa)
        p(:mtpc) = p(:mtpc) + alfa*large_component_values(:mtpc)
        q(:mtpc) = q(:mtpc) + alfa*small_component_values(:mtpc)
        !
        if (method(orbital_index)==2 .or. method(orbital_index)==4) then
            !
            !   Set up right-hand side for variational equations; integrate
            !   variational equations outwards and inwards; store small
            !   component at join point each time
            !
            origin_slope_parameter = 0.0_real64
            call build_variation_source_term (orbital_index)
            radial_equation_type = 3
            call start (orbital_index, radial_equation_type, origin_slope_parameter, pv, small_component_origin_coefficient, qv)
            call integrate_radial_outward (orbital_index, matching_grid_index, pv, qv)
            qjpov = qv(matching_grid_index)
            call integrate_radial_inward (orbital_index, matching_grid_index, pv, qv, mtpv)
            qjpiv = qv(matching_grid_index)
            !
            !   Generate continuous solutions
            !
            mtpvc = max(mtpc,mtpv)
            alfa = -(qjpiv - qjpov)/(qjpih - qjpoh)
            pv(:mtpvc) = pv(:mtpvc) + alfa*large_component_values(:mtpvc)
            qv(:mtpvc) = qv(:mtpvc) + alfa*small_component_values(:mtpvc)
            !
            ta(1) = 0.0_real64
            ta(2:mtpc) = radial_grid_derivative(2:mtpc)*(p(2:mtpc)**2+q(2:mtpc)**2)
            mtp = mtpc
            call quad (dnorm)
            !
            mtp = min(mtpc,mtpvc)
            ta(1) = 0.0_real64
            ta(2:mtp) = radial_grid_derivative(2:mtp)*(p(2:mtp)*pv(2:mtp)+q(2:mtp)*qv(2:mtp))
            call quad (crnorm)
            !
            ta(1) = 0.0_real64
            ta(2:mtpvc) = radial_grid_derivative(2:mtpvc)*(pv(2:mtpvc)**2+qv(2:mtpvc)**2)
            mtp = mtpvc
            call quad (dvnorm)
            !
            !   Determine deleps required to normalize new solution to
            !   first order: modified form of solution to a quadratic
            !   equation (see Press et al.)
            !
            aa = dvnorm
            bb = crnorm + crnorm
            cc = dnorm - 1.0_real64
            discr = bb*bb - 4.0_real64*aa*cc
            if (discr > 0.0_real64) then
                qq = -0.5_real64*(bb + sign(1.0_real64,bb)*sqrt(discr))
                root1 = cc/qq
                root2 = qq/aa
                pmx = 0.0_real64
                apm = 0.0_real64
                do i = 2, matching_grid_index
                    pati = p(i)
                    if (pati > pmx) then
                        pmx = pati
                        loc1 = i
                    end if
                    abpi = abs(pati)
                    if (abpi <= apm) cycle
                    apm = abpi
                    loc2 = i
                end do
                if (pmx /= 0.0_real64) then
                    ratio = apm/abs(pmx)
                    if (ratio < 10.0_real64) then
                        loc = loc1
                    else
                        loc = loc2
                    end if
                else
                    loc = loc2
                end if
                test1 = p(loc) + root1*pv(loc)
                test2 = p(loc) + root2*pv(loc)
                if (test1>0.0_real64 .and. test2<0.0_real64) then
                    dele = root1
                else if (test1<0.0_real64 .and. test2>0.0_real64) then
                    dele = root2
                else if (test1>0.0_real64 .and. test2>0.0_real64) then
                    if (test1 < test2) then
                        dele = root1
                    else
                        dele = root2
                    end if
                else if (test1<0.0_real64 .and. test2<0.0_real64) then
                    if (test1 > test2) then
                        dele = root1
                    else
                        dele = root2
                    end if
                end if
            else
                dele = -bb/(aa + aa)
            end if
            !
            !   Generate new solution
            !
            mtp0 = max(mtpc,mtpvc)
            p0 = p0 + dele*alfa*p0h
            p(2:mtp0) = p(2:mtp0) + dele*pv(2:mtp0)
            q(2:mtp0) = q(2:mtp0) + dele*qv(2:mtp0)
        else
            mtp0 = mtpc
        end if
        !
        !   Debug printout
        !
        if (ldbpr(23)) call print_radial_wavefunction (orbital_index)
        !
        !   Count nodes in large component function; determine sign at first
        !   oscillation, effective quantum number; note that node counting
        !   is never enforced on the small component
        !
        call count (p, mtp0, nnp, first_oscillation_sign)
        !
        !   DEBUG PRINTOUT
        !
        if (ldbpr(22)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(" Iteration number: ",1i2,", limit: ",1i2,/," Present estimate of P0; ",1d21.14,/&
                &," Present estimate of E(J): ",1d21.14,", DELEPS: ",1d21.14,/," Lower bound on energy: ",1d21.14,&
                &", upper bound: ",1d21.14,/," Join point: ",1i4,", Maximum tabulation point:",1i4,/,&
                &" Number of nodes counted: ",1i2,", Correct number: ",1i2,/," Sign of P at first oscillation: ",f3.0)') &
                kount, nloops, p0, orbital_energy_parameters(orbital_index), energy_correction, epsmin, &
                epsmax, matching_grid_index, mtp, nnp, nnodep(orbital_index), first_oscillation_sign
        end if
        !
        !   Proceed according to method
        !
        if (method(orbital_index) > 2) then
            if (check .and. first_oscillation_sign<0.0_real64) then
                orbital_sign_reversal_count = 1
                p0 = -p0
                p(2:mtp0) = -p(2:mtp0)
                q(2:mtp0) = -q(2:mtp0)
            end if
        else
            node_count_difference = nnp - nnodep(orbital_index)
            effective_principal_n = nnp + orbital_l(orbital_index) + 1
            call update_trial_orbital_energy (orbital_index, first_oscillation_sign, effective_principal_n, node_count_difference, &
                & energy_correction, &
                & solve_failed, orbital_sign_reversal_count)
            if (solve_failed) cycle
        end if
        exit
    end do
    !
    !   Solution found
    !
    return
    !
    return
    !
end subroutine solve_radial_orbital
