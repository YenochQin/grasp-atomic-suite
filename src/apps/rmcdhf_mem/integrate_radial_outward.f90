!***********************************************************************
!                                                                      *
subroutine integrate_radial_outward(orbital_index, matching_grid_index, large_component_values, small_component_values)
    ! Module dependencies
    !                                                                      *
    !   This subroutine carries out the step-by-step outward integration   *
    !   of a pair of inhomogeneous Dirac radial equations.                 *
    !                                                                      *
    !   arguments:                                                         *
    !                                                                      *
    !      J:   (Input) Orbital index of function to be computed           *
    !      JP:  (Input) The join point; the outward integration stops      *
    !           at this tabulation index                                   *
    !      P,Q: (Input and output) on input, elements 1 to 3 of both       *
    !           arrays must be tabulated; on output, the arrays are        *
    !           tabulated up to point JP                                   *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford   Last updated: 08 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_radial_points
    use radial_grid_state,           only: transformed_grid_step, radial_derivative_over_radius
    use radial_dirac_workspace,            only: tf ,tg ,xu, xv
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
    !-----------------------------------------------
    !
    !   One global initialization
    !
    dkhb2 = 0.5_real64*transformed_grid_step*dble(orbital_kappa(orbital_index))
    !
    !   Tabulate P(r) and Q(r) by step-by-step integration
    !
    !   Initializations: set quantities for I = 3
    !
    i = 3
    dkhb2f = dkhb2*radial_derivative_over_radius(i)
    cpi = 1.0_real64 + dkhb2f
    cmi = 1.0_real64 - dkhb2f
    pi = large_component_values(i)
    qi = small_component_values(i)
    tfi = tf(i)
    tgi = tg(i)
    !
    !   March out to from I = 4 to I = JP
    !
    !XHH Use doo-loop
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
        tfi = tf(i)
        tgi = tg(i)
        ucim1 = cmim1*pim1 - tfim1*qim1 + xu(i-1)
        udim1 = cpim1*qim1 - tgim1*pim1 + xv(i-1)
        uei = cpi*cmi - tfi*tgi
        pi = (cmi*ucim1 - tfi*udim1)/uei
        qi = (cpi*udim1 - tgi*ucim1)/uei
        large_component_values(i) = pi
        small_component_values(i) = qi
    end do
    !
    return
end subroutine integrate_radial_outward
