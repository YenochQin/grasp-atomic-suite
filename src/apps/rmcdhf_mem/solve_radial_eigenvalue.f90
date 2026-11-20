module solve_radial_eigenvalue_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: solve_radial_eigenvalue
contains

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
        use atomic_calculation_state, only: speed_of_light_au
        use radial_grid_state, only: radial_grid_derivative, radial_derivative_over_radius, transformed_grid_step
        use orbital_csf_basis, only: orbital_kappa
        use scf_radial_potentials, only: large_exchange_source, small_exchange_source, orbital_y_potential
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a, radial_workspace_b
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use radial_quadrature, only: integrate_radial_grid
        use radial_wavefunction_derivatives, only: differentiate_radial_wavefunctions
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
        integration_last_index = orbital_last_grid_index(orbital_index)
        !
        !   Exchange term
        !
        radial_workspace_a(1) = 0.0_real64
        do i = 2, integration_last_index
            radial_workspace_a(i) = (radial_large_components(i,orbital_index)*small_exchange_source(i)-radial_small_components(i, &
                & orbital_index)*large_exchange_source(i))*radial_derivative_over_radius(i)
        end do
        call integrate_radial_grid (piece1)
        piece1 = speed_of_light_au*piece1
        !
        !   Direct term
        !
        radial_workspace_a(1) = 0.0_real64
        do i = 2, integration_last_index
            radial_workspace_a(i) = (radial_large_components(i,orbital_index)**2 + radial_small_components(i, &
                & orbital_index)**2)*orbital_y_potential(i)*radial_derivative_over_radius(i)
        end do
        call integrate_radial_grid (piece2)
        !
        !   Kinetic energy terms
        !
        radial_workspace_a(1) = 0.0_real64
        do i = 2, integration_last_index
            radial_workspace_a(i) = radial_small_components(i,orbital_index)**2*radial_grid_derivative(i)
        end do
        call integrate_radial_grid (piece3)
        piece3 = 2.0_real64*speed_of_light_au*speed_of_light_au*piece3
        !
        radial_workspace_a(1) = 0.0_real64
        do i = 2, integration_last_index
            radial_workspace_a(i) = (radial_large_components(i,orbital_index)*radial_small_components(i, &
                & orbital_index))*radial_derivative_over_radius(i)
        end do
        call integrate_radial_grid (piece4)
        piece4 = -2.0_real64*dble(orbital_kappa(orbital_index))*speed_of_light_au*piece4
        !
        call differentiate_radial_wavefunctions (orbital_index)
        radial_workspace_a(1) = 0.0_real64
        do i = 2, integration_last_index
            radial_workspace_a(i) = radial_large_components(i,orbital_index)*radial_workspace_b(i) - &
                radial_small_components(i,orbital_index)*radial_workspace_a(i)
        end do
        call integrate_radial_grid (piece5)
        piece5 = speed_of_light_au*piece5/transformed_grid_step
        !
        !   Assembly
        !
        solve_radial_eigenvalue = piece1 + piece2 + piece3 + piece4 + piece5
        !
        return
        !
    end function solve_radial_eigenvalue
end module solve_radial_eigenvalue_interface
