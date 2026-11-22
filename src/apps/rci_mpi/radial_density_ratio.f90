! Evaluate the radial density ratio used in the field-shift calculation.
module radial_density_ratio_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: radial_density_ratio
contains

    !    This subprogram computes the overlap of the orbital tabulated in
    !    the arrays  P  and  Q  with maximum tabulation point  MTPO  with
    !    a hydrogenic orbital with parameters  NP  KAPPA  Z .
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 09 Oct 1992
    !
    ! ...Switches:
    real(real64) function radial_density_ratio (large_component_values, small_component_values, last_radial_index, &
        & principal_n, kappa, &
        & effective_nuclear_charge)
        !  Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use radial_grid_state, only: radial_grid_derivative
        use hydrogenic_radial_components, only: hydrogenic_large_component, hydrogenic_small_component
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use dirac_coulomb_radial_orbital, only: evaluate_dirac_coulomb_orbital
        use radial_quadrature, only: integrate_radial_grid
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
        ! Set up the hydrogenic orbital
        call evaluate_dirac_coulomb_orbital (principal_n, kappa, effective_nuclear_charge, eh, pzh, &
            hydrogenic_large_component, hydrogenic_small_component, mtph)
        ! Compute the overlap
        integration_last_index = min(mtph,last_radial_index)
        do i = 2, integration_last_index
            if (radial_grid_derivative(i) > 0.0219_real64) cycle
            k = i
        end do
        integration_last_index = k
        radial_workspace_a(1) = 0.0_real64
        radial_workspace_a(2:integration_last_index) = &
            (large_component_values(2:integration_last_index)*large_component_values(2:integration_last_index)+ &
            small_component_values(2:integration_last_index)* &
            & small_component_values(2:integration_last_index))*radial_grid_derivative(2:integration_last_index)
        ! TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
        call integrate_radial_grid (result)
        radial_workspace_a(2:integration_last_index) = (hydrogenic_large_component(2:integration_last_index)* &
            hydrogenic_large_component(2:integration_last_index)+ &
            hydrogenic_small_component(2:integration_last_index)*hydrogenic_small_component(2:integration_last_index))* &
                radial_grid_derivative(2:integration_last_index)
        call integrate_radial_grid (result1)
        radial_density_ratio = result/result1
        return
    end function radial_density_ratio
end module radial_density_ratio_interface
