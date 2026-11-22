! Solve the homogeneous radial Dirac equation for the initial-orbital estimate.
module solve_homogeneous_radial_equation_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: solve_homogeneous_radial_equation
contains

    !   This routine solves the homogeneous Dirac radial equation.
    !
    !   Arguments:  IORB : (Input) Index of orbital
    !               FAIL : (Output) .TRUE. if solution not obtained
    !
    !   The direct potential is assumed tabulated in the COMMON array YP
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 09 Dec 1992
    subroutine solve_homogeneous_radial_equation(orbital_index, solve_failed)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use atomic_calculation_state, only: working_tolerance, nuclear_charge, speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid_derivative
        use orbital_csf_basis, only: orbital_energy_parameters, orbital_kappa, orbital_principal_n
        use scf_radial_potentials, only: orbital_y_potential
        use radial_quadrature_workspace, only: radial_workspace_a, integration_last_index
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use radial_dirac_coefficients, only: calculate_radial_dirac_coefficients
        use radial_origin_solution, only: solve_radial_origin

        ! Procedure interfaces
        use integrate_homogeneous_radial_segment_interface, only: integrate_homogeneous_radial_segment
        use integrate_radial_tail_interface, only: integrate_radial_tail

        ! Module dependencies
        use radial_node_count, only: count_radial_nodes
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        logical , intent(out) :: solve_failed

        ! Local declarations
        integer :: max_matching_iterations, npiorb, nkiorb, nakabs, nliorb, nnp, nrel, kount, jp, &
            radial_equation_type, nstrt, i, npc
        real(real64) :: epslon, csq, alpha, fkabs, fkap2, zalpha, gamma, ebym, &
            emin, emax, orbital_energy_correction, trial_orbital_energy_parameter, small_component_origin_coefficient, pfjpo, &
                & qfjpo, pfjpi, qfjpi, ratio, first_oscillation_sign, &
            dnorm, small_component_matching_residual, dnfac

        ! Initial values
        data max_matching_iterations/ 75/
        ! A solution is deemed continuous when the relative mismatch in the
        ! small component is within EPSLON
        epslon = working_tolerance*0.1e00_real64
        ! Establish the number of nodes in the large component
        npiorb = orbital_principal_n(orbital_index)
        nkiorb = orbital_kappa(orbital_index)
        nakabs = abs(nkiorb)
        if (nkiorb < 0) then
            nliorb = nakabs - 1
        else
            nliorb = nakabs
        end if
        nnp = npiorb - nliorb - 1
        ! Establish the bounds on, and an estimate of, the eigenvalue
        csq = speed_of_light_au*speed_of_light_au
        alpha = 1.0_real64/speed_of_light_au
        nrel = npiorb - nakabs
        fkabs = dble(nakabs)
        fkap2 = fkabs*fkabs
        zalpha = orbital_y_potential(radial_point_count)*alpha
        if (zalpha < fkabs) then
            gamma = sqrt(fkap2 - zalpha*zalpha)
            ebym = 1.0_real64/sqrt(1.0_real64 + (zalpha/(gamma + nrel + 0.5_real64))**2)
            emin = (1.0_real64 - ebym)*csq
        else
            emin = 0.25_real64*csq/dble(npiorb*npiorb)
        end if
        zalpha = nuclear_charge*alpha
        if (zalpha < fkabs) then
            gamma = sqrt(fkap2 - zalpha*zalpha)
            ebym = 1.0_real64/sqrt(1.0_real64 + (zalpha/(gamma + nrel))**2)
            orbital_energy_parameters(orbital_index) = (1.0_real64 - ebym)*csq
        else
            orbital_energy_parameters(orbital_index) = csq
        end if
        if (zalpha < fkabs) then
            gamma = sqrt(fkap2 - zalpha*zalpha)
            ebym = 1.0_real64/sqrt(1.0_real64 + (zalpha/(gamma + nrel - 0.5_real64))**2)
            emax = (1.0_real64 - ebym)*csq
        else
            emax = csq + csq
        end if
        orbital_energy_correction = 0.0_real64
        ! Initialize
        solve_failed = .false.
        kount = -1
        ! Iteration loop begins here
        do
            kount = kount + 1
            if (kount > max_matching_iterations) then
                solve_failed = .true.
                return
            end if
            ! Generate estimate of eigenvalue for this iteration
            trial_orbital_energy_parameter = orbital_energy_parameters(orbital_index) + orbital_energy_correction
            if (trial_orbital_energy_parameter>emin .and. trial_orbital_energy_parameter<emax) then
                orbital_energy_parameters(orbital_index) = trial_orbital_energy_parameter
            else
                orbital_energy_parameters(orbital_index) = 0.5_real64*(emin + emax)
            end if
            ! E is the positive binding energy. An invalid trial energy cannot
            ! produce the decaying bound-state tail used by TAIL.
            if (.not.ieee_is_finite(orbital_energy_parameters(orbital_index)) .or. &
                orbital_energy_parameters(orbital_index)<=0.0_real64 .or. orbital_energy_parameters(orbital_index) &
                & >=2.0_real64*csq) then
                solve_failed = .true.
                return
            end if
            ! Set up arrays TF and TG; find join point
            call calculate_radial_dirac_coefficients (orbital_index, jp)
            ! Initialize outward integration
            radial_equation_type = 1
            call solve_radial_origin (orbital_index, radial_equation_type, orbital_origin_parameters(orbital_index), &
                & radial_large_components(:, &
                & orbital_index), small_component_origin_coefficient, &
                & radial_small_components(:,orbital_index))
            ! Continue outward integration to classical turning point
            nstrt = 6
            call integrate_homogeneous_radial_segment (orbital_index, nstrt, jp, radial_large_components(:,orbital_index), &
                & radial_small_components(:,orbital_index))
            pfjpo = radial_large_components(jp,orbital_index)
            qfjpo = radial_small_components(jp,orbital_index)
            ! Initialize inward integration
            call integrate_radial_tail (orbital_index, radial_large_components(:,orbital_index), radial_small_components(:, &
                & orbital_index), jp, integration_last_index)
            ! Continue inward integration to classical turning point
            nstrt = integration_last_index
            call integrate_homogeneous_radial_segment (orbital_index, nstrt, jp, radial_large_components(:,orbital_index), &
                & radial_small_components(:,orbital_index))
            pfjpi = radial_large_components(jp,orbital_index)
            qfjpi = radial_small_components(jp,orbital_index)
            ! Make large component continuous, determine mismatch in small
            ! component
            ratio = pfjpo/pfjpi
            radial_large_components(jp:radial_point_count,orbital_index) = radial_large_components(jp:radial_point_count, &
                & orbital_index)*ratio
            radial_small_components(jp:radial_point_count,orbital_index) = radial_small_components(jp:radial_point_count, &
                & orbital_index)*ratio
            ! Count nodes
            call count_radial_nodes (radial_large_components(:,orbital_index), integration_last_index, npc, first_oscillation_sign)
            ! Correct the energy estimate if the number of nodes is wrong
            if (npc > nnp) then
                emin = orbital_energy_parameters(orbital_index)
                orbital_energy_correction = 0.5_real64*(emax - emin)
                cycle
            else if (npc < nnp) then
                emax = orbital_energy_parameters(orbital_index)
                orbital_energy_correction = -0.5_real64*(emax - emin)
                cycle
            end if
            ! Correct number of nodes
            !
            ! Compute 'norm' of solution
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:radial_point_count) = (radial_large_components(2:radial_point_count,orbital_index)**2+ &
                & radial_small_components(2:radial_point_count,orbital_index)**2)*radial_grid_derivative(2:radial_point_count)
            integration_last_index = radial_point_count
            call integrate_radial_grid (dnorm)
            if (.not.ieee_is_finite(dnorm) .or. dnorm<=0.0_real64) then
                solve_failed = .true.
                return
            end if
            ! Determine correction to eigenvalue from magic formula
            ! correct slope at origin
            qfjpi = qfjpi*ratio
            small_component_matching_residual = qfjpi - qfjpo
            if (.not.ieee_is_finite(small_component_matching_residual) .or. &
                .not.ieee_is_finite(qfjpo)) then
                solve_failed = .true.
                return
            end if
            if (abs(small_component_matching_residual) > epslon*abs(qfjpo)) then
                orbital_energy_correction = speed_of_light_au*radial_large_components(jp,orbital_index)* &
                    & small_component_matching_residual/dnorm
                if (.not.ieee_is_finite(orbital_energy_correction)) then
                    solve_failed = .true.
                    return
                end if
                ! A large correction (notably for the first h/i orbitals) can make
                ! the old multiplicative upper bound negative. Keep the current
                ! trial inside the updated interval; the next iteration bisects it
                ! if the proposed correction falls outside the bounds.
                if (orbital_energy_correction < 0.0_real64) then
                    emax = max(orbital_energy_parameters(orbital_index), emax*(1.0_real64 - 0.2e00_real64* &
                        & abs(orbital_energy_correction/ &
                        & orbital_energy_parameters(orbital_index))))
                else
                    emin = min(orbital_energy_parameters(orbital_index), emin*(1.0_real64 + 0.2e00_real64* &
                        & abs(orbital_energy_correction/ &
                        & orbital_energy_parameters(orbital_index))))
                end if
                orbital_origin_parameters(orbital_index) = orbital_origin_parameters(orbital_index)/sqrt(dnorm)
                cycle
            end if
            exit
        end do
        ! Normalize
        dnfac = 1.0_real64/sqrt(dnorm)
        orbital_origin_parameters(orbital_index) = orbital_origin_parameters(orbital_index)*dnfac
        radial_large_components(:radial_point_count,orbital_index) = radial_large_components(:radial_point_count, &
            & orbital_index)*dnfac
        radial_small_components(:radial_point_count,orbital_index) = radial_small_components(:radial_point_count, &
            & orbital_index)*dnfac
        ! Find maximum tabulation point
        i = radial_point_count + 1
        do
            i = i - 1
            if (abs(radial_large_components(i,orbital_index)) < epslon) then
                radial_large_components(i,orbital_index) = 0.0_real64
                radial_small_components(i,orbital_index) = 0.0_real64
            else
                orbital_last_grid_index(orbital_index) = i
                exit
            end if
        end do
        return
    end subroutine solve_homogeneous_radial_equation
end module solve_homogeneous_radial_equation_interface
