! Propagate the homogeneous radial Dirac solution along a segment of the grid.
module integrate_homogeneous_radial_segment_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: integrate_homogeneous_radial_segment
contains

    !   This  subroutine continues the solution of the homogeneous Dirac
    !   radial equation  from tabulation point NSTRT to tabulation point
    !   NEND. The algorithm of J E Sienkiewicz and W E Baylis, J Phys B:
    !   At Mol Phys 20 (1987) 5145-5156, p 5155, is used.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 08 Dec 1992
    subroutine integrate_homogeneous_radial_segment(orbital_index, first_radial_index, last_radial_index, &
        & large_component_values, small_component_values)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: max_radial_points
        use radial_grid_state, only: radial_derivative_over_radius, transformed_grid_step
        use radial_dirac_workspace, only: large_equation_coefficients, small_equation_coefficients
        use orbital_csf_basis, only: orbital_kappa
        use radial_stepping_coefficients, only: starting_step_weights
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: orbital_index
        integer , intent(in) :: first_radial_index
        integer , intent(in) :: last_radial_index
        real(real64) , intent(inout) :: large_component_values(max_radial_points)
        real(real64) , intent(inout) :: small_component_values(max_radial_points)

        ! Local constants
        ! Weight the current and three preceding P and Q samples in the direction of integration.
        ! The starting_step_weights imported above weight the derivatives.
        real(real64), parameter :: current_value_weight = -0.50_real64
        real(real64), parameter :: one_step_back_weight = 0.50_real64
        real(real64), parameter :: two_steps_back_weight = 0.75_real64
        real(real64), parameter :: three_steps_back_weight = 0.25_real64

        ! Local declarations
        integer :: idiff, loc, j, jp1, jm1
        real(real64) :: tbh, tc1, fk, cc, pj, qj, fac, ppj, qpj, pjm1, qjm1, &
            ppjm1, qpjm1, pjm2, qjm2, ppjm2, qpjm2, pjm3, qjm3, ppjm3, qpjm3, pjm4&
            , qjm4, ppjm4, qpjm4, rppj, rqpj, ccrpor, cpjp1, cmjp1, fjp1, gjp1, &
            denom, factor, pjp1, qjp1, ppjp1, qpjp1, pjp2, qjp2, ppjp2, qpjp2, &
            pjp3, qjp3, ppjp3, qpjp3, pjp4, qjp4, ppjp4, qpjp4, rpmj, rqmj, cpjm1&
            , cmjm1, fjm1, gjm1
        real(real64) :: c1, c2, c3, c4, c5, c6
        logical :: out
        ! Initializations
        c1 = starting_step_weights(1)
        c2 = starting_step_weights(2)
        c3 = starting_step_weights(3)
        c4 = starting_step_weights(4)
        c5 = starting_step_weights(5)
        c6 = starting_step_weights(6)
        tbh = 2.0_real64/transformed_grid_step
        tc1 = 2.0_real64*c1
        ! Determine whether integration is inward or outwarD
        idiff = last_radial_index - first_radial_index
        if (idiff > 0) then
            out = .true.
        else if (idiff == 0) then
            return
        else if (idiff < 0) then
            out = .false.
        end if
        ! Overall initializations
        fk = dble(orbital_kappa(orbital_index))
        cc = c1*fk*transformed_grid_step
        ! Perform integration depending on case
        if (out) then
            ! Initialization for outward integration
            loc = first_radial_index
            pj = large_component_values(loc)
            qj = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppj = (-fac*pj) - tbh*large_equation_coefficients(loc)*qj
            qpj = fac*qj - tbh*small_equation_coefficients(loc)*pj
            loc = loc - 1
            pjm1 = large_component_values(loc)
            qjm1 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjm1 = (-fac*pjm1) - tbh*large_equation_coefficients(loc)*qjm1
            qpjm1 = fac*qjm1 - tbh*small_equation_coefficients(loc)*pjm1
            loc = loc - 1
            pjm2 = large_component_values(loc)
            qjm2 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjm2 = (-fac*pjm2) - tbh*large_equation_coefficients(loc)*qjm2
            qpjm2 = fac*qjm2 - tbh*small_equation_coefficients(loc)*pjm2
            loc = loc - 1
            pjm3 = large_component_values(loc)
            qjm3 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjm3 = (-fac*pjm3) - tbh*large_equation_coefficients(loc)*qjm3
            qpjm3 = fac*qjm3 - tbh*small_equation_coefficients(loc)*pjm3
            loc = loc - 1
            pjm4 = large_component_values(loc)
            qjm4 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjm4 = (-fac*pjm4) - tbh*large_equation_coefficients(loc)*qjm4
            qpjm4 = fac*qjm4 - tbh*small_equation_coefficients(loc)*pjm4
            ! March out
            j = first_radial_index - 1
            do
                j = j + 1
                rppj = current_value_weight*pj + one_step_back_weight*pjm1 &
                    + two_steps_back_weight*pjm2 + three_steps_back_weight*pjm3 &
                    + c2*ppj + c3*ppjm1 + c4*ppjm2 + c5*ppjm3 + c6*ppjm4
                rqpj = current_value_weight*qj + one_step_back_weight*qjm1 &
                    + two_steps_back_weight*qjm2 + three_steps_back_weight*qjm3 &
                    + c2*qpj + c3*qpjm1 + c4*qpjm2 + c5*qpjm3 + c6*qpjm4
                jp1 = j + 1
                ccrpor = cc*radial_derivative_over_radius(jp1)
                cpjp1 = 1.0_real64 + ccrpor
                cmjp1 = 1.0_real64 - ccrpor
                fjp1 = tc1*large_equation_coefficients(jp1)
                gjp1 = tc1*small_equation_coefficients(jp1)
                denom = cpjp1*cmjp1 - gjp1*fjp1
                factor = 1.0_real64/denom
                pjp1 = (cmjp1*rppj - fjp1*rqpj)*factor
                qjp1 = (cpjp1*rqpj - gjp1*rppj)*factor
                large_component_values(jp1) = pjp1
                small_component_values(jp1) = qjp1
                if (jp1 >= last_radial_index) exit
                ppjm4 = ppjm3
                qpjm4 = qpjm3
                pjm3 = pjm2
                qjm3 = qjm2
                ppjm3 = ppjm2
                qpjm3 = qpjm2
                pjm2 = pjm1
                qjm2 = qjm1
                ppjm2 = ppjm1
                qpjm2 = qpjm1
                pjm1 = pj
                qjm1 = qj
                ppjm1 = ppj
                qpjm1 = qpj
                pj = pjp1
                qj = qjp1
                fac = fk*radial_derivative_over_radius(jp1)
                ppj = (-fac*pj) - tbh*large_equation_coefficients(jp1)*qj
                qpj = fac*qj - tbh*small_equation_coefficients(jp1)*pj
            end do
        else
            ! Initializations for inward integration
            loc = first_radial_index
            pj = large_component_values(loc)
            qj = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppj = (-fac*pj) - tbh*large_equation_coefficients(loc)*qj
            qpj = fac*qj - tbh*small_equation_coefficients(loc)*pj
            loc = loc + 1
            pjp1 = large_component_values(loc)
            qjp1 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjp1 = (-fac*pjp1) - tbh*large_equation_coefficients(loc)*qjp1
            qpjp1 = fac*qjp1 - tbh*small_equation_coefficients(loc)*pjp1
            loc = loc + 1
            pjp2 = large_component_values(loc)
            qjp2 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjp2 = (-fac*pjp2) - tbh*large_equation_coefficients(loc)*qjp2
            qpjp2 = fac*qjp2 - tbh*small_equation_coefficients(loc)*pjp2
            loc = loc + 1
            pjp3 = large_component_values(loc)
            qjp3 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjp3 = (-fac*pjp3) - tbh*large_equation_coefficients(loc)*qjp3
            qpjp3 = fac*qjp3 - tbh*small_equation_coefficients(loc)*pjp3
            loc = loc + 1
            pjp4 = large_component_values(loc)
            qjp4 = small_component_values(loc)
            fac = fk*radial_derivative_over_radius(loc)
            ppjp4 = (-fac*pjp4) - tbh*large_equation_coefficients(loc)*qjp4
            qpjp4 = fac*qjp4 - tbh*small_equation_coefficients(loc)*pjp4
            ! March in
            j = first_radial_index + 1
            do
                j = j - 1
                rpmj = current_value_weight*pj + one_step_back_weight*pjp1 &
                    + two_steps_back_weight*pjp2 + three_steps_back_weight*pjp3 &
                    - c2*ppj - c3*ppjp1 - c4*ppjp2 - c5*ppjp3 - c6*ppjp4
                rqmj = current_value_weight*qj + one_step_back_weight*qjp1 &
                    + two_steps_back_weight*qjp2 + three_steps_back_weight*qjp3 &
                    - c2*qpj - c3*qpjp1 - c4*qpjp2 - c5*qpjp3 - c6*qpjp4
                jm1 = j - 1
                ccrpor = cc*radial_derivative_over_radius(jm1)
                cpjm1 = 1.0_real64 + ccrpor
                cmjm1 = 1.0_real64 - ccrpor
                fjm1 = tc1*large_equation_coefficients(jm1)
                gjm1 = tc1*small_equation_coefficients(jm1)
                denom = cpjm1*cmjm1 - gjm1*fjm1
                factor = 1.0_real64/denom
                pjm1 = (cpjm1*rpmj + fjm1*rqmj)*factor
                qjm1 = (cmjm1*rqmj + gjm1*rpmj)*factor
                large_component_values(jm1) = pjm1
                small_component_values(jm1) = qjm1
                if (jm1 <= last_radial_index) exit
                ppjp4 = ppjp3
                qpjp4 = qpjp3
                pjp3 = pjp2
                qjp3 = qjp2
                ppjp3 = ppjp2
                qpjp3 = qpjp2
                pjp2 = pjp1
                qjp2 = qjp1
                ppjp2 = ppjp1
                qpjp2 = qpjp1
                pjp1 = pj
                qjp1 = qj
                ppjp1 = ppj
                qpjp1 = qpj
                pj = pjm1
                qj = qjm1
                fac = fk*radial_derivative_over_radius(jm1)
                ppj = (-fac*pj) - tbh*large_equation_coefficients(jm1)*qj
                qpj = fac*qj - tbh*small_equation_coefficients(jm1)*pj
            end do
        end if
        return
    end subroutine integrate_homogeneous_radial_segment
end module integrate_homogeneous_radial_segment_interface
