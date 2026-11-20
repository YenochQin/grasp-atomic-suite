!***********************************************************************
!                                                                      *
subroutine calculate_biorthogonal_overlap
    ! Module dependencies
    !                                                                      *
    !   This subroutine calculates the initial and final state             *
    !   radial overlap matrix                                              *
    !                                                                      *
    !   Written by Per Jonsson                                             *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
        nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
    use radial_wavefunctions, only: mfii, pfii, qfii, mfff, pfff, qfff
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Local declarations
    integer :: i, j, l
    real(real64) :: braket(max_orbitals, max_orbitals)
    real(real64) :: result
    !-----------------------------------------------

    !
    do i = 1, nwii
        do j = 1, nwff
            if (nakii(i) /= nakff(j)) cycle
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

            braket(i,j) = result

            !            WRITE (*,9) '<',NPII(I),NHII(I),'|',NPFF(J),NHFF(J),'> =', &
            !               BRAKET(I,J)
        end do
    end do
    !
    !    9 FORMAT(A,I2,A,A,I2,A,A,E20.13)

    return
end subroutine calculate_biorthogonal_overlap
