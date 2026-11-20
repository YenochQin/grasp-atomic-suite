! Integrate the overlap of two radial spinors weighted by a power of r.
module weighted_radial_overlap
    implicit none (type, external)

    ! Module interface
    private
    public :: rint
contains
    real(real64) function rint (first_orbital, second_orbital, radial_power)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
            radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use radial_quadrature_workspace, only: weighted_integrand => ta, last_radial_index => mtp
        use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
            pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
        use radial_quadrature, only: quad
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital
        integer, intent(in) :: second_orbital
        integer, intent(in) :: radial_power

        ! Local declarations
        integer :: radial_index
        real(real64) :: overlap_integral
        last_radial_index = min(orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        weighted_integrand(1) = 0.0_real64
        do radial_index = 2, last_radial_index
            weighted_integrand(radial_index) = radial_grid(radial_index)**radial_power*(radial_large_components(radial_index, &
                & first_orbital)*radial_large_components(radial_index,second_orbital) + &
                & radial_small_components(radial_index,first_orbital)*radial_small_components(radial_index, &
                & second_orbital))*radial_grid_derivative(radial_index)
        end do
        call quad (overlap_integral)
        rint = overlap_integral
        return
    end function rint
end module weighted_radial_overlap
