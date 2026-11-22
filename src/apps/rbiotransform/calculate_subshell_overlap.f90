! Evaluate radial overlaps for orbitals in one subshell symmetry group.
module calculate_subshell_overlap_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_subshell_overlap
contains

    !   This subroutine calculates the initial and final state
    !   radial overlap matrix
    !
    ! Written by Per Jonsson
    subroutine calculate_subshell_overlap(overlap_matrix, initial_orbital_count, final_orbital_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_grid_derivative
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: final_orbital_last_grid_index, initial_orbital_last_grid_index, &
            final_radial_large_components, final_radial_small_components, initial_radial_large_components, &
            initial_radial_small_components
        use radial_quadrature, only: integrate_radial_grid

        ! Procedure interfaces
        use write_matrix_interface, only: write_matrix
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: initial_orbital_count
        integer, intent(in) :: final_orbital_count
        real(real64), intent(out) :: overlap_matrix(initial_orbital_count,final_orbital_count)

        ! Local declarations
        integer      :: initial_orbital_index, final_orbital_index, radial_index
        real(real64) :: radial_overlap

        do initial_orbital_index = 1, initial_orbital_count
            do final_orbital_index = 1, final_orbital_count
                ! Limit the integral to the supported radial grid points.
                integration_last_index = min(initial_orbital_last_grid_index(initial_orbital_index), &
                    final_orbital_last_grid_index(final_orbital_index))
                ! Tabulate the weighted integrand for radial quadrature. The first grid value is unused.
                radial_workspace_a(1) = 0.0_real64
                do radial_index = 2, integration_last_index
                    radial_workspace_a(radial_index) = &
                        & (initial_radial_large_components(radial_index,initial_orbital_index)* &
                            final_radial_large_components(radial_index,final_orbital_index) + &
                        & initial_radial_small_components(radial_index,initial_orbital_index)* &
                            final_radial_small_components(radial_index,final_orbital_index)) &
                        & *radial_grid_derivative(radial_index)
                end do
                ! Integrate the tabulated radial integrand.
                call integrate_radial_grid (radial_overlap)
                overlap_matrix(initial_orbital_index,final_orbital_index) = radial_overlap
            end do
        end do
        ! Print out
        write (*, *) "********************"
        write (*, *) " S matrix from calculate_subshell_overlap"
        write (*, *) "********************"

        call write_matrix (overlap_matrix, initial_orbital_count, final_orbital_count, initial_orbital_count, final_orbital_count)
        return
    end subroutine calculate_subshell_overlap
end module calculate_subshell_overlap_interface
