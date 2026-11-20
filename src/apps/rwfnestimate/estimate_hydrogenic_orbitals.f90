module estimate_hydrogenic_orbitals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: estimate_hydrogenic_orbitals
contains

    !***********************************************************************
    !                                                                      *
    subroutine estimate_hydrogenic_orbitals(matched_orbital_indices, matched_subshell_count, modify, nradial)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: input_unit
        !                                                                      *
        !   This  subroutine  is  used  to  produce  estimates  of  the wave   *
        !   functions from the hydrogenic approximation. The screening cons-   *
        !   tant is taken to be 0.4 (Changed, see comments below,XHH)          *
        !                                                                      *
        !   Call(s) to: [ERWF]: DCBSRW.                                        *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: nuclear_charge
        use initial_orbital_selection, only: max_orbitals, orbital_initialized
        use orbital_csf_basis, only: orbital_energy_parameters, orbital_principal_n, orbital_kappa, orbital_symmetry_labels
        use radial_wavefunctions, only: radial_large_components, radial_small_components, orbital_origin_parameters, &
            orbital_last_grid_index
        use initial_orbital_source, only: orbital_source_description
        use hydrogenic_screening, only: orbital_screening_constants
        use standard_io_units, only: stderr_unit
        use dirac_coulomb_radial_orbital, only: evaluate_dirac_coulomb_orbital
        implicit none (type, external)

        ! Arguments and local declarations
        integer , intent(in) :: matched_subshell_count
        logical , intent(in) :: modify
        integer , intent(in) :: matched_orbital_indices(max_orbitals)
        integer :: matched_subshell_index, orbital_index
        integer, intent(in) :: nradial
        real(real64) :: effective_nuclear_charge, delz
        !-----------------------------------------------
        !
        !
        delz = 0.0_real64
        if (nradial == 4)  then
            write (stderr_unit, *) "Enter increase in Z for correlation orbitals"
            read  (input_unit, *) delz
        end if
        write (stderr_unit, *) "Orbital Z_eff for hydrogenic orbitals"
        do matched_subshell_index = 1, matched_subshell_count
            orbital_index = matched_orbital_indices(matched_subshell_index)
            if (orbital_initialized(orbital_index)) cycle
            if (modify) then
                write (stderr_unit, *) "Input new value >"
                read (input_unit, *) orbital_screening_constants(orbital_index)
            end if
            effective_nuclear_charge = nuclear_charge - orbital_screening_constants(orbital_index) + delz
            write (stderr_unit, "(I2,A2,F10.2)") orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                & effective_nuclear_charge
            !           ...Calculate radial wavefunctions
            call evaluate_dirac_coulomb_orbital (orbital_principal_n(orbital_index), orbital_kappa(orbital_index), &
                & effective_nuclear_charge, &
                & orbital_energy_parameters(orbital_index), &
                & orbital_origin_parameters(orbital_index), radial_large_components(:,orbital_index), &
                radial_small_components(:,orbital_index), orbital_last_grid_index(orbital_index))
            orbital_initialized(orbital_index) = .true.
            orbital_source_description(orbital_index) = "Hyd"
        end do
        write (stderr_unit, *)
        return
    end subroutine estimate_hydrogenic_orbitals
end module estimate_hydrogenic_orbitals_interface
