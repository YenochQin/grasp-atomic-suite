! Integrate the coupled radial Dirac equations from the origin to the matching point.
module integrate_radial_outward_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: integrate_radial_outward
contains

    ! The large and small components must already be tabulated at the first three points.
    ! Propagate both components to matching_grid_index to join the inward solution.
    ! Written by Farid A Parpia, at Oxford. Last updated: 08 Dec 1992
    subroutine integrate_radial_outward(orbital_index, matching_grid_index, large_component_values, small_component_values)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_radial_points
        use radial_grid_state,           only: transformed_grid_step, radial_derivative_over_radius
        use radial_dirac_workspace,            only: large_equation_coefficients, small_equation_coefficients, &
            large_radial_source, small_radial_source
        use orbital_csf_basis, only: orbital_kappa
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: orbital_index
        integer , intent(in) :: matching_grid_index
        real(real64) , intent(inout) :: large_component_values(max_radial_points)
        real(real64) , intent(inout) :: small_component_values(max_radial_points)

        ! Local declarations
        integer :: i
        real(real64) :: dkhb2, dkhb2f, cpi, cmi, pi, qi, tfi, tgi, cpim1, cmim1, &
            pim1, qim1, tfim1, tgim1, ucim1, udim1, uei
        ! One global initialization
        dkhb2 = 0.5_real64*transformed_grid_step*dble(orbital_kappa(orbital_index))
        ! Tabulate P(r) and Q(r) by step-by-step integration
        !
        ! Initializations: set quantities for I = 3
        i = 3
        dkhb2f = dkhb2*radial_derivative_over_radius(i)
        cpi = 1.0_real64 + dkhb2f
        cmi = 1.0_real64 - dkhb2f
        pi = large_component_values(i)
        qi = small_component_values(i)
        tfi = large_equation_coefficients(i)
        tgi = small_equation_coefficients(i)
        !    March out to from I = 4 to I = JP
        !
        ! XHH Use doo-loop
        do i = 4, matching_grid_index
            cpim1 = cpi
            cmim1 = cmi
            dkhb2f = dkhb2*radial_derivative_over_radius(i)
            cpi = 1.0_real64 + dkhb2f
            cmi = 1.0_real64 - dkhb2f
            pim1 = pi
            qim1 = qi
            tfim1 = tfi
            tgim1 = tgi
            tfi = large_equation_coefficients(i)
            tgi = small_equation_coefficients(i)
            ucim1 = cmim1*pim1 - tfim1*qim1 + large_radial_source(i-1)
            udim1 = cpim1*qim1 - tgim1*pim1 + small_radial_source(i-1)
            uei = cpi*cmi - tfi*tgi
            pi = (cmi*ucim1 - tfi*udim1)/uei
            qi = (cpi*udim1 - tgi*ucim1)/uei
            large_component_values(i) = pi
            small_component_values(i) = qi
        end do
        return
    end subroutine integrate_radial_outward
end module integrate_radial_outward_interface
