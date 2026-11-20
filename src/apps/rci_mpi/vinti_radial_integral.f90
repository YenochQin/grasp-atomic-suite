module vinti_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: vinti_radial_integral
contains

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
        use debug_output_settings, only: rmcdhf_debug_flags
        use radial_grid_state, only: radial_derivative_over_radius, transformed_grid_step
        use orbital_csf_basis, only: orbital_kappa, orbital_principal_n, orbital_symmetry_labels
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a, radial_workspace_b
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use radial_wavefunction_derivatives, only: differentiate_radial_wavefunctions
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital
        integer, intent(in) :: second_orbital

        ! Local declarations
        integer :: i, kpj, kpk, ifact1, ifact2
        real(real64) :: piece1, fact1, fact2, piece2
        !-----------------------------------------------
        !
        integration_last_index = max(orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        !
        !   Piece involving derivatives
        !
        call differentiate_radial_wavefunctions (second_orbital)
        radial_workspace_a(1) = 0.0_real64
        do i = 2, integration_last_index
            radial_workspace_a(i) = radial_large_components(i,first_orbital)*radial_workspace_a(i) + &
                radial_small_components(i,first_orbital)*radial_workspace_b(i)
        end do
        call integrate_radial_grid (piece1)
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
        radial_workspace_a(1) = 0.0_real64
        do i = 2, integration_last_index
            radial_workspace_a(i) = radial_derivative_over_radius(i)*(fact1*radial_large_components(i, &
                first_orbital)*radial_large_components(i, &
                & second_orbital) + &
                & fact2*radial_small_components(i,first_orbital)*radial_small_components(i,second_orbital))
        end do
        call integrate_radial_grid (piece2)
        !
        vinti_radial_integral = piece1 - piece2
        !
        !   Debug printout
        !
        if (rmcdhf_debug_flags(6)) then
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
end module vinti_radial_integral_interface
