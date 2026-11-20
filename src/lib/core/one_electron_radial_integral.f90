! Compute the one-electron radial Dirac matrix element for two orbitals with the same kappa.
module one_electron_radial_integral
    implicit none (type, external)

    ! Module interface
    private
    public :: rinti
contains
    real(real64) function rinti (bra_orbital, ket_orbital, potential_mode)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
        use atomic_calculation_state,           only: speed_of_light_au
        use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
            radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use nuclear_potential_state,          only: negative_radius_times_potential
        use orbital_csf_basis, only: orbital_kappa, orbital_principal_n, orbital_symmetry_labels
        use radial_quadrature_workspace,          only: &
            integrand_or_large_derivative => ta, small_derivative_workspace => tb, last_radial_index => mtp
        use radial_wavefunctions,          only: orbital_last_grid_index,radial_large_components,radial_small_components
        use radial_wavefunction_derivatives, only: dpbdt
        use radial_quadrature, only: quad
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_orbital
        integer, intent(in) :: ket_orbital
        integer, intent(in) :: potential_mode

        ! Local declarations
        integer :: radial_index
        real(real64) :: derivative_term, small_component_mass_term, kappa_term, nuclear_potential_term
        if (orbital_kappa(bra_orbital) /= orbital_kappa(ket_orbital)) then
            write (*, '("RINTI: Attempt to calculate I(",1i2,1a2,",",1i2,1a2,")")') &
                orbital_principal_n(bra_orbital), orbital_symmetry_labels(bra_orbital), &
                & orbital_principal_n(ket_orbital), orbital_symmetry_labels(ket_orbital)
            stop
        end if
        last_radial_index = max(orbital_last_grid_index(bra_orbital),orbital_last_grid_index(ket_orbital))
        ! The derivative workspace contains h*dP/dt and h*dQ/dt. divide the first integral by h.
        call dpbdt (ket_orbital)
        integrand_or_large_derivative(1) = 0.0_real64
        do radial_index = 2, last_radial_index
            integrand_or_large_derivative(radial_index) = radial_small_components(radial_index, &
                & bra_orbital)*integrand_or_large_derivative(radial_index) - radial_large_components(radial_index, &
                & bra_orbital)*small_derivative_workspace(radial_index)
        end do
        call quad (derivative_term)
        derivative_term = speed_of_light_au*derivative_term/transformed_grid_step
        integrand_or_large_derivative(1) = 0.0_real64
        do radial_index = 2, last_radial_index
            integrand_or_large_derivative(radial_index) = &
                & radial_grid_derivative(radial_index)*radial_small_components(radial_index, &
                & bra_orbital)*radial_small_components(radial_index,ket_orbital)
        end do
        call quad (small_component_mass_term)
        small_component_mass_term = -2.0_real64*speed_of_light_au*speed_of_light_au*small_component_mass_term
        integrand_or_large_derivative(1) = 0.0_real64
        do radial_index = 2, last_radial_index
            integrand_or_large_derivative(radial_index) = &
                & radial_derivative_over_radius(radial_index)*(radial_large_components(radial_index, &
                & bra_orbital)*radial_small_components(radial_index,ket_orbital) + &
                & radial_small_components(radial_index,bra_orbital)*radial_large_components(radial_index,ket_orbital))
        end do
        call quad (kappa_term)
        kappa_term = kappa_term*speed_of_light_au*dble(orbital_kappa(ket_orbital))
        ! Only mode zero includes the nuclear potential.
        if (potential_mode == 0) then
            integrand_or_large_derivative(1) = 0.0_real64
            do radial_index = 2, last_radial_index
                integrand_or_large_derivative(radial_index) = &
                    & radial_derivative_over_radius(radial_index)*negative_radius_times_potential(radial_index)* &
                    & (radial_large_components(radial_index,bra_orbital)* &
                    & radial_large_components(radial_index,ket_orbital) + radial_small_components(radial_index, &
                    & bra_orbital)*radial_small_components(radial_index,ket_orbital))
            end do
            call quad (nuclear_potential_term)
            nuclear_potential_term = -nuclear_potential_term
        else
            nuclear_potential_term = 0.0_real64
        end if
        rinti = derivative_term + small_component_mass_term + kappa_term + nuclear_potential_term
        if (potential_mode==0 .and. ldbpr(4)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/," I (",1i2,1a2,",",1i2,1a2,") = ",1p,d19.12,/)') &
                orbital_principal_n(bra_orbital), &
                & orbital_symmetry_labels(bra_orbital), &
                & orbital_principal_n(ket_orbital), orbital_symmetry_labels(ket_orbital), &
                rinti
        end if
        if (potential_mode/=0 .and. ldbpr(5)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/," K (",1i2,1a2,",",1i2,1a2,") = ",1p,d19.12,/)') &
                orbital_principal_n(bra_orbital), &
                & orbital_symmetry_labels(bra_orbital), &
                & orbital_principal_n(ket_orbital), orbital_symmetry_labels(ket_orbital), &
                rinti
        end if
        return
        return
    end function rinti
end module one_electron_radial_integral
