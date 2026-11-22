! Add the off-diagonal integral contributions to the current matrix element.
module apply_off_diagonal_integrals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: apply_off_diagonal_integrals
contains

    !   This  routine  includes  the  contribution from the off-diagonal
    !   I(a,b) integrals in the 'exchange' term.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 08 Dec 1992
    subroutine apply_off_diagonal_integrals
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: speed_of_light_au
        use radial_grid_state, only: radial_grid, radial_derivative_over_radius, transformed_grid_step
        use nuclear_potential_state, only: negative_radius_times_potential
        use orbital_csf_basis, only: orbital_kappa
        use scf_radial_potentials, only: large_exchange_source, small_exchange_source
        use scf_iteration_state, only: direct_coefficient_count, direct_coefficients, direct_orbital_indices
        use radial_quadrature_workspace, only: radial_workspace_a, radial_workspace_b
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use radial_wavefunction_derivatives, only: differentiate_radial_wavefunctions
        implicit none (type, external)

        ! Local declarations
        integer :: k, orbital_index, mfi, i
        real(real64) :: twoc, coeff, fk, rporii, pfi, qfi, zbci

        twoc = speed_of_light_au + speed_of_light_au
        do k = 1, direct_coefficient_count
            orbital_index = direct_orbital_indices(k)
            call differentiate_radial_wavefunctions (orbital_index)
            mfi = orbital_last_grid_index(orbital_index)
            coeff = direct_coefficients(k)
            fk = dble(orbital_kappa(orbital_index))
            do i = 2, mfi
                rporii = 1.0_real64/(transformed_grid_step*radial_derivative_over_radius(i))
                pfi = radial_large_components(i,orbital_index)
                qfi = radial_small_components(i,orbital_index)
                zbci = negative_radius_times_potential(i)/speed_of_light_au
                large_exchange_source(i) = large_exchange_source(i) + coeff*(radial_workspace_a(i)*rporii+ &
                    fk*pfi-(twoc*radial_grid(i)+zbci)*qfi)
                small_exchange_source(i) = small_exchange_source(i) + coeff*(radial_workspace_b(i)*rporii-fk*qfi+zbci*pfi)
            end do
        end do
        return
    end subroutine apply_off_diagonal_integrals
end module apply_off_diagonal_integrals_interface
