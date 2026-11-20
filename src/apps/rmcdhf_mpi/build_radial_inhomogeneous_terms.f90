module build_radial_inhomogeneous_terms_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_radial_inhomogeneous_terms
contains

    !***********************************************************************
    !                                                                      *
    subroutine build_radial_inhomogeneous_terms(orbital_index)
        ! Module dependencies
        !                                                                      *
        !   This  SUBROUTINE  sets  up the arrays XU and XV, for use by  the   *
        !   subprograms IN and  OUT.                                           *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: large_deferred_correction, small_deferred_correction
        use radial_grid_state, only: radial_point_count, radial_derivative_over_radius, transformed_grid_step
        use radial_dirac_workspace, only: large_radial_source, small_radial_source
        use scf_radial_potentials, only: large_exchange_source, small_exchange_source
        implicit none (type, external)

        ! Arguments
        integer  :: orbital_index

        ! Local declarations
        integer :: nm1
        real(real64) :: dmhh
        !-----------------------------------------------
        !
        !   Define constants
        !
        dmhh = -transformed_grid_step*0.5_real64
        !
        !   Set up arrays XU and XV; since XU(1), XV(1) are never used,
        !   set them to some arbitrary value
        !
        nm1 = radial_point_count - 1
        large_radial_source(1) = 0.0_real64
        small_radial_source(1) = 0.0_real64
        large_radial_source(2:nm1) = dmhh*(large_exchange_source(3:nm1+1)*radial_derivative_over_radius(3:nm1+1)+ &
            large_exchange_source(2:nm1)*radial_derivative_over_radius(2:nm1)) + large_deferred_correction(2&
            :nm1)
        small_radial_source(2:nm1) = dmhh*(small_exchange_source(3:nm1+1)*radial_derivative_over_radius(3:nm1+1)+ &
            small_exchange_source(2:nm1)*radial_derivative_over_radius(2:nm1)) + small_deferred_correction(2&
            :nm1)
        !
        large_radial_source(radial_point_count) = &
            dmhh*large_exchange_source(radial_point_count)*radial_derivative_over_radius(radial_point_count) + &
            & large_deferred_correction(radial_point_count)
        small_radial_source(radial_point_count) = &
            dmhh*large_exchange_source(radial_point_count)*radial_derivative_over_radius(radial_point_count) + &
            & small_deferred_correction(radial_point_count)
        !
        return
    end subroutine build_radial_inhomogeneous_terms
end module build_radial_inhomogeneous_terms_interface
