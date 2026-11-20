!***********************************************************************
!                                                                      *
subroutine build_breit_bessel_kernels(orbital_a, orbital_b, ik, iw, k)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This routine evaluates the functions                               *
    !                                                                      *
    !                  (2K+1)!!                                            *
    !      BESSJ  =    --------  J  (w   r) - 1  =  PHI  (w   r) - 1       *
    !                         K   K   ab               K   ab              *
    !                  (w   r)                                             *
    !                    ab                                                *
    !   and                                                                *
    !                         K+1                                          *
    !                  (w   r)                                             *
    !                    ab                                                *
    !      BESSN  =  - --------  N  (w   r) - 1  =  PSI  (w   r) - 1       *
    !                  (2K-1)!!   K   ab               K   ab              *
    !                                                                      *
    !                                                                      *
    !   where J and N ARE spherical Bessel functions, and                  *
    !                                                                      *
    !                w   = ABS ( E(IA) - E(IB) )/c                         *
    !                 ab                                                   *
    !                                                                      *
    !   where  E(I)  is the eigenvalue for orbital  I . The writeup (B J   *
    !   McKenzie, I P Grant, and P H Norrington, Computer Phys Commun 21   *
    !   (1980) 233-246) is incorrect in its description of the output of   *
    !   this routine.                                                      *
    !                                                                      *
    !   The  routine uses equations given in M Abramowitz and I A STegun   *
    !   to evaluate the functions. Devices are used to reduce the number   *
    !   of actual evaluations of these functions.                          *
    !                                                                      *
    !                                           Last update: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: orbital_index_radix
    use bessel_function_workspace, only: bessj, bessn, bj, max_radial_points, tc, td, wij
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use atomic_calculation_state, only: accy, max_radial_points, speed_of_light_au
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use orbital_csf_basis, only: orbital_energy_parameters, orbital_principal_n, orbital_symmetry_labels
    use breit_bessel_cache_keys, only: keep
    use breit_frequency_scale, only: wfact
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_a
    integer, intent(in) :: orbital_b
    integer, intent(in) :: ik
    integer, intent(in) :: iw
    integer, intent(in) :: k

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    integer :: icode, iwkp, ikkp, l, i, nn, radial_index, jchan, irem, iswap
    real(real64) :: epsi, w, wa, xbess1, xbess2, s1, s2, dfnm, dfn, ssn, scn&
        , sn, cn, obwa, b, skeep
    !-----------------------------------------------
    !
    epsi = sqrt(0.1_real64*accy)
    !
    !   Form unique label symmetric in IA, IB
    !
    icode = max(orbital_a,orbital_b) + orbital_label_radix*(min(orbital_a,orbital_b) + orbital_label_radix*k)
    !
    !   Function in position; return
    !
    if (icode == keep(ik,iw)) return
    !
    !   Function not in position; is it available in BESS arrays?
    !
    w = wfact*abs(orbital_energy_parameters(orbital_a)-orbital_energy_parameters(orbital_b))/speed_of_light_au
    wij(iw) = w
    !
    do iwkp = 1, 2
        do ikkp = 1, 2
            if (keep(ikkp,iwkp) /= icode) cycle
            !
            !   Function found move into position
            !
            keep(ik,iw) = icode
            if (ldbpr(7)) then
                call ensure_debug_output()
                write (angular_debug_unit, '(93x,i2,a2,2x,i2,a2,2x,i2,2x,"(",i2,",",i2,")",2x,"(",i2,",",i2,")")') &
                    orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), orbital_principal_n(orbital_b) &
                    & , orbital_symmetry_labels(orbital_b), k, &
                    ikkp, iwkp, ik, iw
            end if
            bessj(ik,iw,:radial_point_count) = bessj(ikkp,iwkp,:radial_point_count)
            bessn(ik,iw,:radial_point_count) = bessn(ikkp,iwkp,:radial_point_count)
            return
        end do
    end do
    !
    !   Function not found; evaluate it
    !
    if (ldbpr(7)) then
        call ensure_debug_output()
        write (angular_debug_unit, '(93x,i2,a2,2x,i2,a2,2x,i2,2x,"New",6x,"(",i2,",",i2,")")') &
            orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), orbital_principal_n(orbital_b), &
            & orbital_symmetry_labels(orbital_b), k, ik, iw
    end if
    !
    keep(ik,iw) = icode
    !
    if (w < epsi**2) then
        !
        !   Negligible w
        !
        bessj(ik,iw,:radial_point_count) = 0.0_real64
        bessn(ik,iw,:radial_point_count) = 0.0_real64
        return
        !
    end if
    !
    nn = k
    !
    bessj(ik,iw,1) = 0.0_real64
    bessn(ik,iw,1) = 0.0_real64
    !
    !   Use a four-term power series for low w*r
    !
    do radial_index = 2, radial_point_count
        wa = -0.5_real64*(radial_grid(radial_index)*w)**2
        xbess1 = 1.0_real64
        xbess2 = 1.0_real64
        s1 = 0.0_real64
        s2 = 0.0_real64
        do i = 1, 4
            xbess1 = xbess1*wa/dble(i*(2*(nn + i) + 1))
            xbess2 = xbess2*wa/dble(i*(2*(i - nn) - 1))
            s1 = s1 + xbess1
            s2 = s2 + xbess2
            if (abs(xbess1)>=abs(s1)*epsi .or. abs(xbess2)>=abs(s2)*epsi) &
                cycle
            exit
        end do
        if (i > 4) exit
        bessj(ik,iw,radial_index) = s1
        bessn(ik,iw,radial_index) = s2
    end do
    !
    !   If here then calculated whole array using four-term power
    !   series.  Hence return
    !
    if (radial_index > radial_point_count) return
    jchan = radial_index
    !
    !   Use sin/cos expansion when power series requires more than
    !   four terms terms to converge
    !
    if (nn == 0) then
        dfnm = 1.0_real64
        dfn = 1.0_real64
    else
        dfnm = 1.0_real64
        do i = 3, 2*nn - 1, 2
            dfnm = dfnm*dble(i)
        end do
        dfn = dfnm*dble(2*nn + 1)
    end if
    dfnm = 1.0_real64/dfnm
    !
    irem = mod(nn,4)
    !
    select case (irem)
    case (1)
        !
        !   NN = 1, 5, 9, ...
        !
        ssn = -1.0_real64
        scn = 1.0_real64
        iswap = 1
    !
    case (2)
        !
        !   N = 2, 6, 10, ....
        !
        ssn = -1.0_real64
        scn = -1.0_real64
        iswap = 0
    !
    case (3)
        !
        !   NN = 3, 7, 11,...
        !
        ssn = 1.0_real64
        scn = -1.0_real64
        iswap = 1
    !
    case default
        !
        !   NN = 0, 4, 8,...
        !
        ssn = 1.0_real64
        scn = 1.0_real64
        iswap = 0
        !
    end select
    !
    do radial_index = jchan, radial_point_count
        wa = w*radial_grid(radial_index)
        if (iswap == 0) then
            sn = ssn*sin(wa)
            cn = scn*cos(wa)
        else
            sn = ssn*cos(wa)
            cn = scn*sin(wa)
        end if
        obwa = 1.0_real64/wa
        b = obwa
        s1 = b*sn
        s2 = b*cn
        do i = 1, nn
            skeep = sn
            sn = cn
            cn = -skeep
            b = b*obwa*dble((nn + i)*(nn - i + 1))/dble(2*i)
            s1 = s1 + b*sn
            s2 = s2 + b*cn
        end do
        s1 = s1*dfn/wa**nn - 1.0_real64
        s2 = s2*wa**(nn + 1)*dfnm - 1.0_real64
        bessj(ik,iw,radial_index) = s1
        bessn(ik,iw,radial_index) = s2
    end do
    return
    !
    return
    !
end subroutine build_breit_bessel_kernels
