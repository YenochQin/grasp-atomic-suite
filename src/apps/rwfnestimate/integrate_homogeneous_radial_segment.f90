

!***********************************************************************
!                                                                      *
subroutine integrate_homogeneous_radial_segment(orbital_index, first_radial_index, last_radial_index, &
    & large_component_values, small_component_values)
    ! Module dependencies
    !                                                                      *
    !   This  subroutine continues the solution of the homogeneous Dirac   *
    !   radial equation  from tabulation point NSTRT to tabulation point   *
    !   NEND. The algorithm of J E Sienkiewicz and W E Baylis, J Phys B:   *
    !   At Mol Phys 20 (1987) 5145-5156, p 5155, is used.                  *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters, only: max_radial_points, radial_workspace_capacity, max_orbitals
    !USE DEF_C
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_dirac_workspace, only: tf, tg
    use orbital_csf_basis, only: max_orbitals, orbital_kappa
    use radial_stepping_coefficients, only: c
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: orbital_index
    integer , intent(in) :: first_radial_index
    integer , intent(in) :: last_radial_index
    real(real64) , intent(inout) :: large_component_values(max_radial_points)
    real(real64) , intent(inout) :: small_component_values(max_radial_points)

    ! Local constants
    real(real64), parameter :: b1 = -0.50_real64
    real(real64), parameter :: b2 = 0.50_real64
    real(real64), parameter :: b3 = 0.75_real64
    real(real64), parameter :: b4 = 0.25_real64

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
    !-----------------------------------------------
    !
    !
    !   Initializations
    !
    c1 = c(1)
    c2 = c(2)
    c3 = c(3)
    c4 = c(4)
    c5 = c(5)
    c6 = c(6)
    tbh = 2.0_real64/transformed_grid_step
    tc1 = 2.0_real64*c1
    !
    !   Determine whether integration is inward or outwarD
    !
    idiff = last_radial_index - first_radial_index
    if (idiff > 0) then
        out = .true.
    else if (idiff == 0) then
        return
    else if (idiff < 0) then
        out = .false.
    end if
    !
    !   Overall initializations
    !
    fk = dble(orbital_kappa(orbital_index))
    cc = c1*fk*transformed_grid_step
    !
    !   Perform integration depending on case
    !
    if (out) then
        !
        !   Initialization for outward integration
        !
        loc = first_radial_index
        pj = large_component_values(loc)
        qj = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppj = (-fac*pj) - tbh*tf(loc)*qj
        qpj = fac*qj - tbh*tg(loc)*pj
        !
        loc = loc - 1
        pjm1 = large_component_values(loc)
        qjm1 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjm1 = (-fac*pjm1) - tbh*tf(loc)*qjm1
        qpjm1 = fac*qjm1 - tbh*tg(loc)*pjm1
        !
        loc = loc - 1
        pjm2 = large_component_values(loc)
        qjm2 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjm2 = (-fac*pjm2) - tbh*tf(loc)*qjm2
        qpjm2 = fac*qjm2 - tbh*tg(loc)*pjm2
        !
        loc = loc - 1
        pjm3 = large_component_values(loc)
        qjm3 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjm3 = (-fac*pjm3) - tbh*tf(loc)*qjm3
        qpjm3 = fac*qjm3 - tbh*tg(loc)*pjm3
        !
        loc = loc - 1
        pjm4 = large_component_values(loc)
        qjm4 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjm4 = (-fac*pjm4) - tbh*tf(loc)*qjm4
        qpjm4 = fac*qjm4 - tbh*tg(loc)*pjm4
        !
        !   March out
        !
        j = first_radial_index - 1
        do
            j = j + 1
            !
            rppj = b1*pj + b2*pjm1 + b3*pjm2 + b4*pjm3 + c2*ppj + c3*ppjm1 + c4*&
                ppjm2 + c5*ppjm3 + c6*ppjm4
            rqpj = b1*qj + b2*qjm1 + b3*qjm2 + b4*qjm3 + c2*qpj + c3*qpjm1 + c4*&
                qpjm2 + c5*qpjm3 + c6*qpjm4
            !
            jp1 = j + 1
            ccrpor = cc*radial_derivative_over_radius(jp1)
            cpjp1 = 1.0_real64 + ccrpor
            cmjp1 = 1.0_real64 - ccrpor
            fjp1 = tc1*tf(jp1)
            gjp1 = tc1*tg(jp1)
            denom = cpjp1*cmjp1 - gjp1*fjp1
            factor = 1.0_real64/denom
            pjp1 = (cmjp1*rppj - fjp1*rqpj)*factor
            qjp1 = (cpjp1*rqpj - gjp1*rppj)*factor
            large_component_values(jp1) = pjp1
            small_component_values(jp1) = qjp1
            !
            if (jp1 >= last_radial_index) exit
            !
            ppjm4 = ppjm3
            qpjm4 = qpjm3
            !
            pjm3 = pjm2
            qjm3 = qjm2
            ppjm3 = ppjm2
            qpjm3 = qpjm2
            !
            pjm2 = pjm1
            qjm2 = qjm1
            ppjm2 = ppjm1
            qpjm2 = qpjm1
            !
            pjm1 = pj
            qjm1 = qj
            ppjm1 = ppj
            qpjm1 = qpj
            !
            pj = pjp1
            qj = qjp1
            fac = fk*radial_derivative_over_radius(jp1)
            ppj = (-fac*pj) - tbh*tf(jp1)*qj
            qpj = fac*qj - tbh*tg(jp1)*pj
            !
        end do
    else
        !
        !   Initializations for inward integration
        !
        loc = first_radial_index
        pj = large_component_values(loc)
        qj = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppj = (-fac*pj) - tbh*tf(loc)*qj
        qpj = fac*qj - tbh*tg(loc)*pj
        !
        loc = loc + 1
        pjp1 = large_component_values(loc)
        qjp1 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjp1 = (-fac*pjp1) - tbh*tf(loc)*qjp1
        qpjp1 = fac*qjp1 - tbh*tg(loc)*pjp1
        !
        loc = loc + 1
        pjp2 = large_component_values(loc)
        qjp2 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjp2 = (-fac*pjp2) - tbh*tf(loc)*qjp2
        qpjp2 = fac*qjp2 - tbh*tg(loc)*pjp2
        !
        loc = loc + 1
        pjp3 = large_component_values(loc)
        qjp3 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjp3 = (-fac*pjp3) - tbh*tf(loc)*qjp3
        qpjp3 = fac*qjp3 - tbh*tg(loc)*pjp3
        !
        loc = loc + 1
        pjp4 = large_component_values(loc)
        qjp4 = small_component_values(loc)
        fac = fk*radial_derivative_over_radius(loc)
        ppjp4 = (-fac*pjp4) - tbh*tf(loc)*qjp4
        qpjp4 = fac*qjp4 - tbh*tg(loc)*pjp4
        !
        !   March in
        !
        j = first_radial_index + 1
        do
            j = j - 1
            !
            rpmj = b1*pj + b2*pjp1 + b3*pjp2 + b4*pjp3 - c2*ppj - c3*ppjp1 - c4*&
                ppjp2 - c5*ppjp3 - c6*ppjp4
            rqmj = b1*qj + b2*qjp1 + b3*qjp2 + b4*qjp3 - c2*qpj - c3*qpjp1 - c4*&
                qpjp2 - c5*qpjp3 - c6*qpjp4
            !
            jm1 = j - 1
            ccrpor = cc*radial_derivative_over_radius(jm1)
            cpjm1 = 1.0_real64 + ccrpor
            cmjm1 = 1.0_real64 - ccrpor
            fjm1 = tc1*tf(jm1)
            gjm1 = tc1*tg(jm1)
            denom = cpjm1*cmjm1 - gjm1*fjm1
            factor = 1.0_real64/denom
            pjm1 = (cpjm1*rpmj + fjm1*rqmj)*factor
            qjm1 = (cmjm1*rqmj + gjm1*rpmj)*factor
            large_component_values(jm1) = pjm1
            small_component_values(jm1) = qjm1
            !
            if (jm1 <= last_radial_index) exit
            !
            ppjp4 = ppjp3
            qpjp4 = qpjp3
            !
            pjp3 = pjp2
            qjp3 = qjp2
            ppjp3 = ppjp2
            qpjp3 = qpjp2
            !
            pjp2 = pjp1
            qjp2 = qjp1
            ppjp2 = ppjp1
            qpjp2 = qpjp1
            !
            pjp1 = pj
            qjp1 = qj
            ppjp1 = ppj
            qpjp1 = qpj
            !
            pj = pjm1
            qj = qjm1
            fac = fk*radial_derivative_over_radius(jm1)
            ppj = (-fac*pj) - tbh*tf(jm1)*qj
            qpj = fac*qj - tbh*tg(jm1)*pj
            !
        end do
    end if
    !
    return
end subroutine integrate_homogeneous_radial_segment
