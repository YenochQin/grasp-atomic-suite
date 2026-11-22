! Reset the radial source arrays before assembling the next orbital equation.
module clear_radial_inhomogeneous_terms_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: clear_radial_inhomogeneous_terms
contains

    !   This subprogram sets the inhomogeneous terms to zero.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 17 Dec 1992
    subroutine clear_radial_inhomogeneous_terms(orbital_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_point_count
        use radial_dirac_workspace, only: large_radial_source, small_radial_source
        implicit none (type, external)

        ! Arguments
        integer  :: orbital_index

        large_radial_source(:radial_point_count) = 0.0_real64
        small_radial_source(:radial_point_count) = 0.0_real64
        return
    end subroutine clear_radial_inhomogeneous_terms
end module clear_radial_inhomogeneous_terms_interface
