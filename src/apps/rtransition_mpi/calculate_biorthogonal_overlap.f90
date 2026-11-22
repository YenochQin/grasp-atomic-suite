! Form overlaps between the initial and final orbital sets for the biorthogonal transformation.
module calculate_biorthogonal_overlap_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_biorthogonal_overlap
contains

    !   This subroutine calculates the initial and final state
    !   radial overlap matrix
    !
    ! Written by Per Jonsson
    subroutine calculate_biorthogonal_overlap
        ! Module dependencies
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

        do i = 1, initial_orbital_count
            do j = 1, final_orbital_count
                if (initial_orbital_kappa(i) /= final_orbital_kappa(j)) cycle
                ! Limit the integral to the supported radial grid points.
                integration_last_index = min(initial_orbital_last_grid_index(i),final_orbital_last_grid_index(j))
                ! Tabulate the weighted integrand for radial quadrature. The first grid value is unused.
                radial_workspace_a(1) = 0.0_real64
                do l = 2, integration_last_index
                    radial_workspace_a(l) = (initial_radial_large_components(l,i)*final_radial_large_components(l,j) + &
                        initial_radial_small_components(l,i)*final_radial_small_components(l,j))*radial_grid_derivative(l)
                end do
                ! Integrate the tabulated radial integrand.
                call integrate_radial_grid (result)

                braket(i,j) = result

            end do
        end do

        return
    end subroutine calculate_biorthogonal_overlap
end module calculate_biorthogonal_overlap_interface
