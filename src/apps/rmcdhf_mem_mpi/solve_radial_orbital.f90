module solve_radial_orbital_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: solve_radial_orbital
contains

    !***********************************************************************
    !                                                                      *
    subroutine solve_radial_orbital(orbital_index, solve_failed, orbital_sign_reversal_count, matching_grid_index, &
        & radial_node_count)
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
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state, only: speed_of_light_au, minimum_radial_solver_iterations
        use radial_grid_state, only: radial_grid_derivative
        use radial_dirac_workspace, only: trial_large_component, trial_small_component, trial_large_origin_parameter, &
            trial_small_origin_parameter, trial_orbital_last_grid_index
        use orbital_inversion_flags, only: suppress_radial_inversion
        use expected_orbital_nodes, only: expected_radial_node_counts
        use orbital_csf_basis, only: orbital_energy_parameters, orbital_kappa, orbital_l, orbital_principal_n, &
            orbital_symmetry_labels
        use scf_iteration_state, only: radial_solution_methods, lower_trial_energy_bound, upper_trial_energy_bound, &
            initial_lower_trial_energy_bound, initial_upper_trial_energy_bound, asymptotic_nuclear_charge
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: orbital_origin_parameters
        use scf_radial_potentials, only: max_radial_points

        ! Procedure interfaces
        use estimate_orbital_energy_interface, only: estimate_orbital_energy
        use solve_radial_eigenvalue_interface, only: solve_radial_eigenvalue

        ! Module dependencies
        use dirac_coulomb_radial_orbital, only: evaluate_dirac_coulomb_orbital
        use radial_dirac_coefficients, only: calculate_radial_dirac_coefficients

        ! Procedure interfaces
        use clear_radial_inhomogeneous_terms_interface, only: clear_radial_inhomogeneous_terms

        ! Module dependencies
        use radial_origin_solution, only: solve_radial_origin

        ! Procedure interfaces
        use integrate_radial_outward_interface, only: integrate_radial_outward
        use integrate_radial_inward_interface, only: integrate_radial_inward
        use build_radial_inhomogeneous_terms_interface, only: build_radial_inhomogeneous_terms

        ! Module dependencies
        use radial_quadrature, only: integrate_radial_grid

        ! Procedure interfaces
        use build_variation_source_term_interface, only: build_variation_source_term
        use print_radial_wavefunction_interface, only: print_radial_wavefunction

        ! Module dependencies
        use radial_node_count, only: count_radial_nodes

        ! Procedure interfaces
        use update_trial_orbital_energy_interface, only: update_trial_orbital_energy
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        integer, intent(inout) :: orbital_sign_reversal_count
        integer, intent(out) :: matching_grid_index
        integer, intent(out) :: radial_node_count
        logical, intent(out) :: solve_failed

        ! Local declarations
        integer :: solve_iteration, orbital_relativistic_kappa, maximum_solve_iterations, radial_equation_type, &
            & homogeneous_last_grid_index, inhomogeneous_last_grid_index, i, combined_last_grid_index, &
            & variation_last_grid_index, combined_variation_last_grid_index, &
            loc1, loc2, loc, node_count_difference, effective_principal_n
        real(real64) :: large_component_values(max_radial_points), small_component_values(max_radial_points), &
            & pv(max_radial_points), qv(max_radial_points)
        real(real64) :: twocsq, elast, tenej, energy_correction, first_oscillation_sign, homogeneous_outward_match_small, &
            & homogeneous_inward_match_small, inhomogeneous_outward_match_small, &
            inhomogeneous_inward_match_small, dnorm, alfa, p0h, origin_slope_parameter, small_component_origin_coefficient, &
                & qjpov, qjpiv, crnorm, dvnorm, aa, &
            bb, cc, discr, qq, root1, root2, pmx, apm, pati, abpi, ratio, test1, &
            test2, dele
        logical :: check
        !-----------------------------------------------
        !
        !   Initialization
        !
        check = .not.suppress_radial_inversion(orbital_index)
        solve_failed = .false.
        solve_iteration = 0
        orbital_relativistic_kappa = orbital_kappa(orbital_index)
        !     TWOCSQ = 2.0D 00*137.036**2
        twocsq = 2.0_real64*speed_of_light_au*speed_of_light_au
        !
        !   Debug header
        !
        if (rmcdhf_debug_flags(22)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,/," Debug printout active; orbital: ",1i2,1a2)') &
                orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index)
        end if
        !
        maximum_solve_iterations = max(minimum_radial_solver_iterations,3*orbital_principal_n(orbital_index))
        !
        call estimate_orbital_energy (orbital_index)
        elast = orbital_energy_parameters(orbital_index)
        orbital_energy_parameters(orbital_index) = solve_radial_eigenvalue(orbital_index)
        !
        !   Checks on lower bounds
        !
        if (orbital_energy_parameters(orbital_index) < lower_trial_energy_bound) then
            if (orbital_energy_parameters(orbital_index) > 0.0_real64) then
                lower_trial_energy_bound = orbital_energy_parameters(orbital_index)
            else
                write (*, '(" E(",1i2,1a2,") = ",1p,d11.4,"; adjusted to EPSMIN")') &
                    orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), orbital_energy_parameters( &
                    & orbital_index)
                orbital_energy_parameters(orbital_index) = lower_trial_energy_bound
                if (abs(initial_lower_trial_energy_bound - elast)<=1.0e-6_real64 .and. &
                    radial_solution_methods(orbital_index)<=2) then
                    call evaluate_dirac_coulomb_orbital (orbital_principal_n(orbital_index), orbital_relativistic_kappa, &
                        asymptotic_nuclear_charge, &
                        & orbital_energy_parameters(orbital_index), trial_large_origin_parameter, &
                        & trial_large_component, trial_small_component, &
                        & trial_orbital_last_grid_index)
                    write (*, '(" Returned hydrogenic function for ",1i2,1a2," with"," effective charge ",f7.3)') &
                        orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), asymptotic_nuclear_charge
                    return
                end if
            end if
        end if
        !
        !   Check on upper bound
        !
        if (radial_solution_methods(orbital_index) <= 2) then
            if (orbital_energy_parameters(orbital_index) > upper_trial_energy_bound) then
                tenej = 10.0_real64*orbital_energy_parameters(orbital_index)
                upper_trial_energy_bound = min(tenej,twocsq)
                initial_upper_trial_energy_bound = upper_trial_energy_bound
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
            solve_iteration = solve_iteration + 1
            !
            trial_large_origin_parameter = orbital_origin_parameters(orbital_index)
            !
            if (solve_iteration > 1) then
                !
                !   Check that bounds are ordered correctly
                !
                if (upper_trial_energy_bound <= lower_trial_energy_bound) write (*, &
                    '(" Warning: difficulty with node-counting procedure"/,&
                    &" lower bound on energy (",1p,d11.4,") exceeds upper"," bound (",1d11.4,"; E(",1i2,1a2,") = ",1d11.4)') &
                    lower_trial_energy_bound, upper_trial_energy_bound, orbital_principal_n(orbital_index), &
                        orbital_symmetry_labels(orbital_index), &
                    & orbital_energy_parameters(orbital_index&
                    )
                !
                if (solve_iteration>maximum_solve_iterations .or. &
                    & upper_trial_energy_bound-lower_trial_energy_bound<1.0_real64/dble(orbital_principal_n(orbital_index))**3) then
                    write (*, '(" Method ",1i1," unable to solve for ",1i2,1a2," orbital")') &
                        radial_solution_methods(orbital_index), orbital_principal_n(orbital_index), &
                            orbital_symmetry_labels(orbital_index)
                    write (*, '(" Iteration number: ",1i2,", limit: ",1i2,/," Present estimate of P0; ",1d21.14,/,&
                        &" Present estimate of E(J): ",1d21.14,", DELEPS: ",1d21.14,/," Lower bound on energy: ",1d21.14,&
                        &", upper bound: ",1d21.14,/," Join point: ",1i4,", Maximum tabulation point:",1i4,/,&
                        &" Number of nodes counted: ",1i2,", Correct number: ",1i2,/," Sign of P at first oscillation: ",&
                        &f3.0)') solve_iteration - 1, maximum_solve_iterations, trial_large_origin_parameter, &
                            orbital_energy_parameters(orbital_index), &
                            & energy_correction, lower_trial_energy_bound, upper_trial_energy_bound&
                        , matching_grid_index, integration_last_index, radial_node_count, &
                            expected_radial_node_counts(orbital_index), &
                            & first_oscillation_sign
                    solve_failed = .true.
                    return
                end if
            end if
            !
            !   Set up arrays TF and TG; find join point
            !
            call calculate_radial_dirac_coefficients (orbital_index, matching_grid_index)
            !
            !   Set right-hand side to zero to form homogeneous equations;
            !   integrate homogeneous equations outwards and inwards; store
            !   small component at join point each time
            !
            call clear_radial_inhomogeneous_terms (orbital_index)
            radial_equation_type = 1
            call solve_radial_origin (orbital_index, radial_equation_type, trial_large_origin_parameter, &
                large_component_values, trial_small_origin_parameter, small_component_values)
            call integrate_radial_outward (orbital_index, matching_grid_index, large_component_values, small_component_values)
            homogeneous_outward_match_small = small_component_values(matching_grid_index)
            call integrate_radial_inward (orbital_index, matching_grid_index, large_component_values, small_component_values, &
                & homogeneous_last_grid_index)
            homogeneous_inward_match_small = small_component_values(matching_grid_index)
            !
            !   Set up right-hand side for inhomogeneous equations; integrate
            !   inhomogeneous equations outwards and inwards; store small
            !   component at join point each time
            !
            call build_radial_inhomogeneous_terms (orbital_index)
            radial_equation_type = 2
            call solve_radial_origin (orbital_index, radial_equation_type, trial_large_origin_parameter, &
                trial_large_component, trial_small_origin_parameter, trial_small_component)
            call integrate_radial_outward (orbital_index, matching_grid_index, trial_large_component, trial_small_component)
            inhomogeneous_outward_match_small = trial_small_component(matching_grid_index)
            call integrate_radial_inward (orbital_index, matching_grid_index, trial_large_component, trial_small_component, &
                inhomogeneous_last_grid_index)
            inhomogeneous_inward_match_small = trial_small_component(matching_grid_index)
            !
            !   Determine energy adjustment for methods 1 and 2
            !
            if (radial_solution_methods(orbital_index) <= 2) then
                radial_workspace_a(1) = 0.0_real64
                radial_workspace_a(2:inhomogeneous_last_grid_index) = &
                    & (trial_large_component(2:inhomogeneous_last_grid_index)**2 + &
                        trial_small_component(2:inhomogeneous_last_grid_index)**2) &
                    & *radial_grid_derivative(2:inhomogeneous_last_grid_index)
                integration_last_index = inhomogeneous_last_grid_index
                call integrate_radial_grid (dnorm)
                !        DELEPS = 137.036*P(JP)*(QJPII-QJPOI)/DNORM
                energy_correction = speed_of_light_au*trial_large_component(matching_grid_index)* &
                    (inhomogeneous_inward_match_small - &
                    & inhomogeneous_outward_match_small)/dnorm
            end if
            !
            !   Generate the continuous solution
            !
            combined_last_grid_index = max(homogeneous_last_grid_index,inhomogeneous_last_grid_index)
            alfa = -(inhomogeneous_inward_match_small - inhomogeneous_outward_match_small)/(homogeneous_inward_match_small - &
                & homogeneous_outward_match_small)
            p0h = trial_large_origin_parameter
            trial_large_origin_parameter = trial_large_origin_parameter*(1.0_real64 + alfa)
            trial_large_component(:combined_last_grid_index) = trial_large_component(:combined_last_grid_index) + alfa* &
                large_component_values(:combined_last_grid_index)
            trial_small_component(:combined_last_grid_index) = trial_small_component(:combined_last_grid_index) + alfa* &
                small_component_values(:combined_last_grid_index)
            !
            if (radial_solution_methods(orbital_index)==2 .or. radial_solution_methods(orbital_index)==4) then
                !
                !   Set up right-hand side for variational equations; integrate
                !   variational equations outwards and inwards; store small
                !   component at join point each time
                !
                origin_slope_parameter = 0.0_real64
                call build_variation_source_term (orbital_index)
                radial_equation_type = 3
                call solve_radial_origin (orbital_index, radial_equation_type, origin_slope_parameter, pv, &
                    & small_component_origin_coefficient, qv)
                call integrate_radial_outward (orbital_index, matching_grid_index, pv, qv)
                qjpov = qv(matching_grid_index)
                call integrate_radial_inward (orbital_index, matching_grid_index, pv, qv, variation_last_grid_index)
                qjpiv = qv(matching_grid_index)
                !
                !   Generate continuous solutions
                !
                combined_variation_last_grid_index = max(combined_last_grid_index,variation_last_grid_index)
                alfa = -(qjpiv - qjpov)/(homogeneous_inward_match_small - homogeneous_outward_match_small)
                pv(:combined_variation_last_grid_index) = pv(:combined_variation_last_grid_index) + &
                    & alfa*large_component_values(:combined_variation_last_grid_index)
                qv(:combined_variation_last_grid_index) = qv(:combined_variation_last_grid_index) + &
                    & alfa*small_component_values(:combined_variation_last_grid_index)
                !
                radial_workspace_a(1) = 0.0_real64
                radial_workspace_a(2:combined_last_grid_index) = &
                    & radial_grid_derivative(2:combined_last_grid_index) &
                    & *(trial_large_component(2:combined_last_grid_index)**2 + trial_small_component(2:combined_last_grid_index)**2)
                integration_last_index = combined_last_grid_index
                call integrate_radial_grid (dnorm)
                !
                integration_last_index = min(combined_last_grid_index,combined_variation_last_grid_index)
                radial_workspace_a(1) = 0.0_real64
                radial_workspace_a(2:integration_last_index) = &
                    radial_grid_derivative(2:integration_last_index)*(trial_large_component(2:integration_last_index)* &
                        pv(2:integration_last_index)+ &
                    trial_small_component(2:integration_last_index)*qv(2:integration_last_index))
                call integrate_radial_grid (crnorm)
                !
                radial_workspace_a(1) = 0.0_real64
                radial_workspace_a(2:combined_variation_last_grid_index) = &
                    & radial_grid_derivative(2:combined_variation_last_grid_index) &
                    & *(pv(2:combined_variation_last_grid_index)**2 + qv(2:combined_variation_last_grid_index)**2)
                integration_last_index = combined_variation_last_grid_index
                call integrate_radial_grid (dvnorm)
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
                        pati = trial_large_component(i)
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
                    test1 = trial_large_component(loc) + root1*pv(loc)
                    test2 = trial_large_component(loc) + root2*pv(loc)
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
                trial_orbital_last_grid_index = max(combined_last_grid_index,combined_variation_last_grid_index)
                trial_large_origin_parameter = trial_large_origin_parameter + dele*alfa*p0h
                trial_large_component(2:trial_orbital_last_grid_index) = &
                    trial_large_component(2:trial_orbital_last_grid_index) + dele*pv(2:trial_orbital_last_grid_index)
                trial_small_component(2:trial_orbital_last_grid_index) = &
                    trial_small_component(2:trial_orbital_last_grid_index) + dele*qv(2:trial_orbital_last_grid_index)
            else
                trial_orbital_last_grid_index = combined_last_grid_index
            end if
            !
            !   Debug printout
            !
            if (rmcdhf_debug_flags(23)) call print_radial_wavefunction (orbital_index)
            !
            !   Count nodes in large component function; determine sign at first
            !   oscillation, effective quantum number; note that node counting
            !   is never enforced on the small component
            !
            call count_radial_nodes (trial_large_component, trial_orbital_last_grid_index, radial_node_count, &
                & first_oscillation_sign)
            !
            !   DEBUG PRINTOUT
            !
            if (rmcdhf_debug_flags(22)) then
                call ensure_debug_output()
                write (angular_debug_unit, '(" Iteration number: ",1i2,", limit: ",1i2,/," Present estimate of P0; ",1d21.14,/&
                    &," Present estimate of E(J): ",1d21.14,", DELEPS: ",1d21.14,/," Lower bound on energy: ",1d21.14,&
                    &", upper bound: ",1d21.14,/," Join point: ",1i4,", Maximum tabulation point:",1i4,/,&
                    &" Number of nodes counted: ",1i2,", Correct number: ",1i2,/," Sign of P at first oscillation: ",f3.0)') &
                    solve_iteration, maximum_solve_iterations, trial_large_origin_parameter, &
                        orbital_energy_parameters(orbital_index), &
                        & energy_correction, lower_trial_energy_bound, &
                    upper_trial_energy_bound, matching_grid_index, integration_last_index, radial_node_count, &
                        expected_radial_node_counts(orbital_index), &
                        & first_oscillation_sign
            end if
            !
            !   Proceed according to method
            !
            if (radial_solution_methods(orbital_index) > 2) then
                if (check .and. first_oscillation_sign<0.0_real64) then
                    orbital_sign_reversal_count = 1
                    trial_large_origin_parameter = -trial_large_origin_parameter
                    trial_large_component(2:trial_orbital_last_grid_index) = -trial_large_component(2:trial_orbital_last_grid_index)
                    trial_small_component(2:trial_orbital_last_grid_index) = -trial_small_component(2:trial_orbital_last_grid_index)
                end if
            else
                node_count_difference = radial_node_count - expected_radial_node_counts(orbital_index)
                effective_principal_n = radial_node_count + orbital_l(orbital_index) + 1
                call update_trial_orbital_energy (orbital_index, first_oscillation_sign, effective_principal_n, &
                    & node_count_difference, &
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
end module solve_radial_orbital_interface
