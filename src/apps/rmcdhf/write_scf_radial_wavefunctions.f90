module write_scf_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_scf_radial_wavefunctions
contains

    !***********************************************************************
    subroutine write_scf_radial_wavefunctions(output_wavefunction_path)
        ! Module dependencies
        !***********************************************************************
        use radial_grid_state, only: radial_grid
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_kappa, orbital_principal_n
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_output_unit
        character , intent(in) :: output_wavefunction_path*(*)
        integer :: orbital_index, last_radial_index, radial_index
        !-----------------------------------------------
        !
        open(newunit=wavefunction_output_unit, file=output_wavefunction_path, status="UNKNOWN", form="UNFORMATTED", position=&
            "asis")
        write (wavefunction_output_unit) "G92RWF"
        do orbital_index = 1, orbital_count
            last_radial_index = orbital_last_grid_index(orbital_index)
            write (wavefunction_output_unit) orbital_principal_n(orbital_index), orbital_kappa(orbital_index), &
                & orbital_energy_parameters(orbital_index), last_radial_index
            write (wavefunction_output_unit) orbital_origin_parameters(orbital_index), &
                & (radial_large_components(radial_index,orbital_index), &
                & radial_index=1,last_radial_index), &
                (radial_small_components(radial_index,orbital_index),radial_index=1,last_radial_index)
            write (wavefunction_output_unit) &
                & (radial_grid(radial_index),radial_index=1,last_radial_index)               ! This is a waste of resources
        end do
        close(wavefunction_output_unit)

        return
    end subroutine write_scf_radial_wavefunctions
end module write_scf_radial_wavefunctions_interface
