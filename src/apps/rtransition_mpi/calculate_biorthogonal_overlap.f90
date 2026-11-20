module calculate_biorthogonal_overlap_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_biorthogonal_overlap
contains

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
        use radial_grid_state, only: radial_grid_derivative
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use biorthogonal_orbital_sets, only: max_orbitals, final_orbital_count, initial_orbital_count, final_orbital_kappa, &
            initial_orbital_kappa
        use radial_wavefunctions, only: initial_orbital_last_grid_index, initial_radial_large_components, &
            initial_radial_small_components, final_orbital_last_grid_index, final_radial_large_components, &
            final_radial_small_components
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Local declarations
        integer :: i, j, l
        real(real64) :: braket(max_orbitals, max_orbitals)
        real(real64) :: result
        !-----------------------------------------------

        !
        do i = 1, initial_orbital_count
            do j = 1, final_orbital_count
                if (initial_orbital_kappa(i) /= final_orbital_kappa(j)) cycle
                !
                !   Determine the maximum tabulation point for the integrand
                !
                integration_last_index = min(initial_orbital_last_grid_index(i),final_orbital_last_grid_index(j))
                !
                !   Tabulate the integrand as required for SUBROUTINE QUAD; the
                !   value at the first tabulation point is arbitrary
                !
                radial_workspace_a(1) = 0.0_real64
                do l = 2, integration_last_index
                    radial_workspace_a(l) = (initial_radial_large_components(l,i)*final_radial_large_components(l,j) + &
                        initial_radial_small_components(l,i)*final_radial_small_components(l,j))*radial_grid_derivative(l)
                end do
                !
                !   Perform the quadrature
                !
                call integrate_radial_grid (result)

                braket(i,j) = result

                !            WRITE (*,9) '<',NPII(I),NHII(I),'|',NPFF(J),NHFF(J),'> =', &
                !               BRAKET(I,J)
            end do
        end do
        !
        !    9 FORMAT(A,I2,A,A,I2,A,A,E20.13)

        return
    end subroutine calculate_biorthogonal_overlap
end module calculate_biorthogonal_overlap_interface
