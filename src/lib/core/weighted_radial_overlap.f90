! Integrate the overlap of two radial spinors weighted by a power of r.
module weighted_radial_overlap
    implicit none (type, external)

    ! Module interface
    private
    public :: integrate_weighted_radial_overlap
contains
    real(real64) function integrate_weighted_radial_overlap (first_orbital, second_orbital, radial_power)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_grid, radial_grid_derivative
        use radial_quadrature_workspace, only: radial_workspace_a, integration_last_index
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital
        integer, intent(in) :: second_orbital
        integer, intent(in) :: radial_power

        ! Local declarations
        integer :: radial_index
        real(real64) :: overlap_integral
        integration_last_index = min(orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        ! Fill workspace A with the weighted radial overlap integrand.
        radial_workspace_a(1) = 0.0_real64
        do radial_index = 2, integration_last_index
            radial_workspace_a(radial_index) = radial_grid(radial_index)**radial_power*(radial_large_components(radial_index, &
                & first_orbital)*radial_large_components(radial_index,second_orbital) + &
                & radial_small_components(radial_index,first_orbital)*radial_small_components(radial_index, &
                & second_orbital))*radial_grid_derivative(radial_index)
        end do
        call integrate_radial_grid (overlap_integral)
        integrate_weighted_radial_overlap = overlap_integral
        return
    end function integrate_weighted_radial_overlap
end module weighted_radial_overlap
