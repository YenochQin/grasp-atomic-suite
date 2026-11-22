! Tabulate the radial Bessel functions for the selected transition multipole.
module build_transition_bessel_functions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_transition_bessel_functions
contains

    ! This routine evaluates Bessel fuctions J K ( W*R/C ) at the grid
    ! points for K=L-1,L,L+1 and stores  them in the  arrays BJ(..,1),
    ! BJ(..,2),BJ(..,3) respectively. It uses a power series expansion
    ! for  small  r and switches to sin/cos expansion when more than 4
    ! terms in power series are required.
    ! Last revision: 28 Dec 1992
    subroutine build_transition_bessel_functions(w)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use bessel_function_workspace, only: transition_spherical_bessel_values
        use debug_output_settings, only: rmcdhf_debug_flags
        use radial_grid_state, only: radial_point_count, radial_grid
        use transition_matrix_workspace, only: transition_multipole_type, transition_multipole_rank
        use atomic_calculation_state, only: speed_of_light_au, vacuum_speed_of_light_au
        implicit none (type, external)

        ! Arguments
        real(real64) , intent(in) :: w

        ! Local declarations
        integer :: i1, iw, nn, iend, iprod, i, jchan, radial_index, iswap, modnn4, &
            jj, jjj, ii
        real(real64) :: epsi, dfn, wr, wa, xbess1, s1, ssn, scn, sn, cn, b, skeep

        epsi = 1.0e-5_real64
        iw = 1
        nn = transition_multipole_rank - 1
        if (transition_multipole_type /= 0) then
            iw = 2
            nn = transition_multipole_rank
        end if
        do
            iend = 2*nn + 1
            iprod = 1
            do i = 1, iend, 2
                iprod = iprod*i
            end do
            dfn = iprod
            jchan = radial_point_count
            transition_spherical_bessel_values(1,iw) = 1.0_real64
            do radial_index = 2, radial_point_count
                wr = w*radial_grid(radial_index)
                wa = -wr*wr*0.5_real64
                xbess1 = 1.0_real64
                s1 = 1.0_real64
                do i = 1, 4
                    xbess1 = xbess1*wa/dble(i*(2*(nn + i) + 1))
                    s1 = s1 + xbess1
                    if (abs(xbess1) < abs(s1)*epsi) exit
                end do
                if (i > 4) then
                    jchan = radial_index
                    exit
                end if
                transition_spherical_bessel_values(radial_index,iw) = s1*wr**nn/dfn
            end do
            ! Use sin/cos expansion when power series takes longer
            ! than 4 terms to converge
            if (jchan < radial_point_count) then
                iswap = 0
                modnn4 = mod(nn - 1,4) + 1
                select case (modnn4)
                ! NN = 1, 5, 9, ...
                case default
                    ssn = -1.0_real64
                    scn = 1.0_real64
                    iswap = 1
                ! N = 2, 6, 10,....
                case (2)
                    ssn = -1.0_real64
                    scn = -1.0_real64
                ! NN = 3, 7, 11,...
                case (3)
                    ssn = 1.0_real64
                    scn = -1.0_real64
                    iswap = 1
                ! NN = 0, 4, 8,...
                case (4)
                    ssn = 1.0_real64
                    scn = 1.0_real64
                end select
                do radial_index = jchan, radial_point_count
                    wa = w*radial_grid(radial_index)
                    if (iswap <= 0) then
                        sn = ssn*sin(wa)
                        cn = scn*cos(wa)
                    else
                        sn = ssn*cos(wa)
                        cn = scn*sin(wa)
                    end if
                    i = -1
                    s1 = 0.0_real64
                    i = i + 1
                    i1 = i
                    do i = i1, nn
                        if (i == 0) then
                            b = 1.0_real64/wa
                        else
                            b = b*dble((nn + i)*(nn - i + 1))/dble(2*i)/wa
                        end if
                        s1 = s1 + b*sn
                        skeep = sn
                        sn = cn
                        cn = -skeep
                    end do
                    transition_spherical_bessel_values(radial_index,iw) = s1
                end do
            end if
            if (nn >= transition_multipole_rank + 1 .or. transition_multipole_type == 1) exit
            nn = nn + 1
            iw = iw + 1
        end do
        ! Print out Bessel functions if (debug) option set
        if (rmcdhf_debug_flags(16)) then
            do jj = 1, 3
                jjj = transition_multipole_rank - 2 + jj
                call ensure_debug_output()
                write (angular_debug_unit, '(/," Bessel function of order ",i3,/(1p,7d18.10))') &
                    jjj, (transition_spherical_bessel_values(ii,jj),ii=1,radial_point_count)
            end do
        end if
        !   All done
        do jj = 1, 3
            jjj = transition_multipole_rank - 2 + jj
            transition_spherical_bessel_values(:radial_point_count,jj) = &
                transition_spherical_bessel_values(:radial_point_count,jj)*(speed_of_light_au/vacuum_speed_of_light_au)**jjj
        end do

        return
    end subroutine build_transition_bessel_functions
end module build_transition_bessel_functions_interface
