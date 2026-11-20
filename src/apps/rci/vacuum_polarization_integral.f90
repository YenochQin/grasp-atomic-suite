module vacuum_polarization_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: vacuum_polarization_integral
contains

    !***********************************************************************
    !                                                                      *
    real(real64)  function vacuum_polarization_integral (first_orbital,second_orbital)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   Computes nuclear vacuum polarization integrals.                    *
        !                                                                      *
        !   Call(s) to: [LIB92]: QUAD.                                         *
        !                                           Last update: 09 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use debug_output_settings, only: rmcdhf_debug_flags
        use nuclear_charge_distribution, only: zdist
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: first_orbital, second_orbital

        ! Local constants

        ! Local declarations
        integer :: k
        !------------------------------------------------
        !
        integration_last_index = min (orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        !TA(1) = 0.0D 00  !! Why would this be wrong?!!
        radial_workspace_a(1)  = 0
        do k = 2,integration_last_index
            radial_workspace_a(k) = (radial_large_components(k,first_orbital)*radial_large_components(k,second_orbital)+ &
                radial_small_components(k, &
                & first_orbital)* &
                & radial_small_components(k,second_orbital))                  &
                *zdist(k)
        end do
        call integrate_radial_grid (vacuum_polarization_integral)
        !
        if (rmcdhf_debug_flags(9)) then
            call ensure_debug_output()
            write (angular_debug_unit,'(/"VPINTF: V (",1i2,1a2,",",1i2,1a2,") = ",1pd19.12)') &
                orbital_principal_n(first_orbital),orbital_symmetry_labels(first_orbital),orbital_principal_n(second_orbital), &
                & orbital_symmetry_labels(second_orbital),vacuum_polarization_integral
        end if
        !
        return
        !
        !
    end function vacuum_polarization_integral
end module vacuum_polarization_integral_interface
