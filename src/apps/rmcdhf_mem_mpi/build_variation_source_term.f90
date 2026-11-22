! Build the source term used to update the radial solution when its energy changes.
module build_variation_source_term_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_variation_source_term
contains

    !   This  subprogram  sets up the inhomogeneous terms for the varia-
    !   tion equations.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 17 Dec 1992
    subroutine build_variation_source_term(orbital_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid_derivative, transformed_grid_step
        use radial_dirac_workspace, only: large_radial_source, small_radial_source
        use radial_wavefunctions, only: radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments
        integer  :: orbital_index

        ! Local declarations
        integer :: radial_index
        real(real64) :: hhc

        hhc = 0.5_real64*transformed_grid_step/speed_of_light_au
        ! Set up arrays TF and TG
        do radial_index = 1, radial_point_count
            large_radial_source(radial_index) = -radial_small_components(radial_index,orbital_index)*hhc* &
                radial_grid_derivative(radial_index)
            small_radial_source(radial_index) = radial_large_components(radial_index,orbital_index)*hhc* &
                radial_grid_derivative(radial_index)
        end do
        large_radial_source(:radial_point_count-1) = large_radial_source(2:radial_point_count) + &
            large_radial_source(:radial_point_count-1)
        small_radial_source(:radial_point_count-1) = small_radial_source(2:radial_point_count) + &
            small_radial_source(:radial_point_count-1)
        return
    end subroutine build_variation_source_term
end module build_variation_source_term_interface
