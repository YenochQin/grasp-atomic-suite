module hydrogenic_radial_overlap_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: hydrogenic_radial_overlap
contains

    !***********************************************************************
    !                                                                      *
    real(real64) function hydrogenic_radial_overlap (large_component_values, small_component_values, last_radial_index, &
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
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use radial_grid_state, only: radial_grid_derivative
        use hydrogenic_radial_components, only: hydrogenic_large_component, hydrogenic_small_component
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use dirac_coulomb_radial_orbital, only: evaluate_dirac_coulomb_orbital
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: last_radial_index
        integer  :: principal_n
        integer  :: kappa
        real(real64)  :: effective_nuclear_charge
        real(real64), intent(in) :: large_component_values(max_radial_points)
        real(real64), intent(in) :: small_component_values(max_radial_points)

        ! Local declarations
        integer :: mtph
        real(real64) :: eh, pzh, result
        !-----------------------------------------------
        !
        !   Set up the hydrogenic orbital
        !
        call evaluate_dirac_coulomb_orbital (principal_n, kappa, effective_nuclear_charge, eh, pzh, &
            hydrogenic_large_component, hydrogenic_small_component, mtph)
        !
        !   Compute the overlap
        !
        integration_last_index = min(mtph,last_radial_index)
        radial_workspace_a(1) = 0.0_real64
        radial_workspace_a(2:integration_last_index) = (large_component_values(2:integration_last_index)* &
            hydrogenic_large_component(2:integration_last_index)+ &
            small_component_values(2:integration_last_index)*hydrogenic_small_component(2:integration_last_index))* &
                radial_grid_derivative(2:integration_last_index)
        call integrate_radial_grid (result)
        !
        hydrogenic_radial_overlap = result
        !
        return
    end function hydrogenic_radial_overlap
end module hydrogenic_radial_overlap_interface
