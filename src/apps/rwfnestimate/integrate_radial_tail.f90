

!***********************************************************************
!                                                                      *
subroutine integrate_radial_tail(orbital_index, large_component_values, small_component_values, matching_grid_index, mtp)
    ! Module dependencies
    !                                                                      *
    !   This subroutine begins the inward integration of the homogeneous   *
    !   Dirac radial equation. With only minor modifications, the series   *
    !   given by  J E Sienkiewicz  and W E Baylis, J Phys B: At Mol Phys   *
    !   20 (1987) 5145-5156, p 5155, is used.                              *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters, only: max_radial_points
    use atomic_calculation_state, only: accy, speed_of_light_au, expmax
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use orbital_csf_basis, only: orbital_energy_parameters, orbital_kappa
    use scf_radial_potentials, only: yp
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: orbital_index
    integer , intent(in) :: matching_grid_index
    integer , intent(out) :: mtp
    real(real64) , intent(inout) :: large_component_values(max_radial_points)
    real(real64) , intent(out) :: small_component_values(max_radial_points)

    ! Local declarations
    integer :: radial_index, nm4, m, loc
    real(real64) :: eps, bige, bigebc, bgebc2, beta, fk, rjp, qem, offset, t&
        , zloc, gam2, zlbb, dnu, fkmzbb, sump, sumq, em, cm, emfact, ovltrm, &
        pterm, qterm, exptrm
    !-----------------------------------------------
    !
    !
    !   Initialize
    !
    eps = accy*0.1_real64
    bige = -orbital_energy_parameters(orbital_index)
    bigebc = bige/speed_of_light_au
    bgebc2 = bigebc/speed_of_light_au
    beta = sqrt((-bige*(2.0_real64 + bgebc2)))
    fk = dble(orbital_kappa(orbital_index))
    !
    !   Find MTP
    !
    radial_index = matching_grid_index
    nm4 = radial_point_count - 4
    rjp = radial_grid(matching_grid_index)
    qem = 0.25_real64*expmax
    do
        radial_index = radial_index + 1
        if (radial_index <= nm4) then
            if (.not.(beta*(radial_grid(radial_index)-rjp) > qem)) cycle
            mtp = radial_index
        else
            write (*, '("TAIL: Grid may be of insufficient extent")')
            mtp = nm4
        end if
        exit
    end do
    !
    !   Compute offset for exponential function
    !
    offset = beta*radial_grid(mtp)
    !
    !   Tabulate tail points
    !
    loc = radial_point_count + 1
    do radial_index = mtp, radial_point_count
        !
        t = 2.0_real64*beta*radial_grid(radial_index)
        !
        zloc = yp(radial_index)
        gam2 = fk**2 - (zloc/speed_of_light_au)**2
        zlbb = zloc/beta
        dnu = zlbb*(1.0_real64 + bgebc2)
        fkmzbb = fk - zlbb
        !
        m = -1
        sump = 0.0_real64
        sumq = 0.0_real64
        m = m + 1
        em = dble(m)
        if (m == 0) then
            cm = 1.0_real64
            emfact = 1.0_real64
        else
            cm = cm*(gam2 - (dnu - em)**2)
            emfact = emfact*em
        end if
        ovltrm = cm*(t**(dnu - em)/emfact)
        pterm = ovltrm*(fkmzbb + em)*beta
        qterm = ovltrm*(fkmzbb - em)*bigebc
        sump = sump + pterm
        sumq = sumq + qterm
        do while(abs(pterm/sump)>=eps .or. abs(qterm/sumq)>=eps)
            m = m + 1
            em = dble(m)
            if (m == 0) then
                cm = 1.0_real64
                emfact = 1.0_real64
            else
                cm = cm*(gam2 - (dnu - em)**2)
                emfact = emfact*em
            end if
            ovltrm = cm*(t**(dnu - em)/emfact)
            pterm = ovltrm*(fkmzbb + em)*beta
            qterm = ovltrm*(fkmzbb - em)*bigebc
            sump = sump + pterm
            sumq = sumq + qterm
        end do
        exptrm = exp((-0.5_real64*t) + offset)
        large_component_values(radial_index) = sump*exptrm
        small_component_values(radial_index) = sumq*exptrm
        if (large_component_values(radial_index) /= 0.0_real64) cycle
        loc = radial_index + 1
        exit
    end do
    large_component_values(loc:radial_point_count) = 0.0_real64
    small_component_values(loc:radial_point_count) = 0.0_real64
    !
    return
    !
    return
    !
end subroutine integrate_radial_tail
