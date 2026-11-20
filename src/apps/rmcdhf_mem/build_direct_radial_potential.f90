!***********************************************************************
!                                                                      *
subroutine build_direct_radial_potential(orbital_index)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This subroutine  tabulates the  potential function Y(r) (Eq (14)   *
    !   in  I P Grant, B J McKenzie, P H Norrington, D F Mayers, and N C   *
    !   Pyper, Computer  Phys  Commun  21 (1980) 211) for orbital J. The   *
    !   function is tabulated in the COMMON array YP.                      *
    !                                                                      *
    !   Call(s) to: [LIB92]: DRAW, YZK                                     *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
    !   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: orbital_index_radix
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_potential_state, only: first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
    use orbital_csf_basis, only: max_orbitals, orbital_principal_n, orbital_symmetry_labels
    use scf_radial_potentials, only: max_radial_points, xp, xq, yp
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use hartree_yz_functions, only: yzk
    use radial_ascii_plot, only: draw
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_index

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    integer :: i, index, label, k, ioy1, ioy2, nb3, row_count, ii, ii1, ii2, ii3
    real(real64) :: coeff
    !-----------------------------------------------
    !   Debug printout: composition

    if (ldbpr(29) .or. ldbpr(30)) then
        call ensure_debug_output()
        write (angular_debug_unit, '(/,/," Direct potential for ",1i2,1a2," orbital :"/,/)') &
            orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index)
    end if
    !
    !   Initialize array YP with the nuclear potential piece
    !
    !   Since YA() below  contains contributions from THIS node only,
    !   the initialization should be in consistence with that.

    yp(:radial_point_count) = negative_radius_times_potential(:radial_point_count)

    do index = 1, nycof

        ! Decode information in label
        label = nya(index)
        k = mod(label,orbital_label_radix)
        label = label/orbital_label_radix
        ioy1 = mod(label,orbital_label_radix)
        ioy2 = label/orbital_label_radix
        coeff = ya(index)

        if (ldbpr(29)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,25x,"(",1i2,")"/,1x,1p,d21.14,"* Y    (",1i2,1a2,",",1i2,1a2,")")') &
                k, coeff, orbital_principal_n(ioy1), orbital_symmetry_labels(ioy1), &
                & orbital_principal_n(ioy2)&
                , orbital_symmetry_labels(ioy2)
        end if

        call yzk (k, ioy1, ioy2)                ! Accumulate contributions
        yp(:radial_point_count) = yp(:radial_point_count) - coeff*tb(:radial_point_count)
    end do
    !
    !   Debug printout
    !
    if (ldbpr(30)) then
        call ensure_debug_output()
        write (angular_debug_unit, '(/,/,3(" --------- r --------- ------- Y (r) -------"))')
        nb3 = radial_point_count/3
        if (3*nb3 == radial_point_count) then
            row_count = nb3
        else
            row_count = nb3 + 1
        end if
        do ii = 1, row_count
            ii1 = ii
            ii2 = ii1 + row_count
            ii3 = ii2 + row_count
            if (ii3 <= radial_point_count) then
                call ensure_debug_output()
                write (angular_debug_unit, '(1p,6(1x,1d21.14))') &
                    radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2), radial_grid(ii3), yp(ii3&
                    )
            else if (ii2 <= radial_point_count) then
                call ensure_debug_output()
                write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2)
            else
                call ensure_debug_output()
                write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), yp(ii1)
            end if
        end do
        call draw (yp, 1.0_real64, yp, 0.0_real64, radial_point_count)
    end if
    !
    return
    !
    return
    !
end subroutine build_direct_radial_potential
