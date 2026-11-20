!***********************************************************************
!                                                                      *
subroutine calculate_subshell_overlap(s, nshli, nshlf)
    ! Module dependencies
    !                                                                      *
    !   This subroutine calculates the initial and final state             *
    !   radial overlap matrix                                              *
    !                                                                      *
    !   Written by Per Jonsson                                             *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use orbital_csf_basis, only: max_orbitals
    use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
        nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use radial_quadrature, only: quad

    ! Procedure interfaces
    use write_matrix_interface, only: write_matrix
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: nshli
    integer, intent(in) :: nshlf
    real(real64)  :: s(nshli,nshlf)

    ! Local declarations
    integer      :: i, j, l
    real(real64) :: result
    !-----------------------------------------------
    !
    do i = 1, nshli
        do j = 1, nshlf
            !
            !   Determine the maximum tabulation point for the integrand
            !
            mtp = min(mfii(i),mfff(j))
            !
            !   Tabulate the integrand as required for SUBROUTINE QUAD; the
            !   value at the first tabulation point is arbitrary
            !
            ta(1) = 0.0_real64
            do l = 2, mtp
                ta(l) = (pfii(l,i)*pfff(l,j) + qfii(l,i)*qfff(l,j))*radial_grid_derivative(l)
            end do
            !
            !   Perform the quadrature
            !
            call quad (result)
            s(i,j) = result
        end do
    end do
    !
    !   Print out
    !
    write (*, *) "********************"
    write (*, *) " S matrix from GETS"
    write (*, *) "********************"

    call write_matrix (s, nshli, nshlf, nshli, nshlf)
    !
    return
end subroutine calculate_subshell_overlap
