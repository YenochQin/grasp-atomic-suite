! Evaluate the double-orbital Slater radial integral using the Hartree Y function.
module slater_radial_integral
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_slater_radial_integral
contains
    real(real64) function calculate_slater_radial_integral (orbital_a, orbital_b, orbital_c, orbital_d, tensor_rank)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use debug_output_settings, only: rmcdhf_debug_flags
        use radial_grid_state, only: radial_derivative_over_radius
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use radial_quadrature_workspace, only: radial_workspace_a, radial_workspace_b, integration_last_index
        use radial_wavefunctions,          only: orbital_last_grid_index, radial_small_components, radial_large_components
        use hartree_yz_functions, only: evaluate_hartree_yz_functions
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_a
        integer, intent(in) :: orbital_b
        integer, intent(in) :: orbital_c
        integer, intent(in) :: orbital_d
        integer, intent(in) :: tensor_rank

        ! Local declarations
        integer :: radial_index
        real(real64) :: slater_integral
        ! Form Y^k from the (b,d) pair in workspace B. Use workspace A for the weighted (a,c) radial integrand.
        call evaluate_hartree_yz_functions (tensor_rank, orbital_b, orbital_d)
        if (tensor_rank==0 .and. orbital_b==orbital_d) then
            integration_last_index = min(orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
        else
            integration_last_index = min(min(integration_last_index,orbital_last_grid_index(orbital_a)), &
                orbital_last_grid_index(orbital_c))
        end if
        radial_workspace_a(1) = 0.0_real64
        do radial_index = 2, integration_last_index
            radial_workspace_a(radial_index) = (radial_large_components(radial_index, &
                & orbital_a)*radial_large_components(radial_index,orbital_c) + radial_small_components(radial_index,orbital_a)* &
                & radial_small_components(radial_index, &
                & orbital_c))*radial_derivative_over_radius(radial_index)*radial_workspace_b(radial_index)
        end do
        call integrate_radial_grid (slater_integral)
        calculate_slater_radial_integral = slater_integral
        if (rmcdhf_debug_flags(10)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"  (",1i1,")",/," R   (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,1a2,") ","= ",&
                &1p,d19.12)') tensor_rank, orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), &
                & orbital_principal_n(orbital_b), &
                & orbital_symmetry_labels(orbital_b), orbital_principal_n(orbital_c)&
                , orbital_symmetry_labels(orbital_c), orbital_principal_n(orbital_d), orbital_symmetry_labels(orbital_d), &
                    & calculate_slater_radial_integral
        end if
        return
        return
    end function calculate_slater_radial_integral
end module slater_radial_integral
