! Save the initial and final biorthogonal radial orbitals for the transition calculation.
module write_transformed_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_transformed_radial_wavefunctions
contains

    !  This subroutine outputs the transformed radial orbitals
    !
    ! Written by Per Jonsson
    subroutine write_transformed_radial_wavefunctions(name)
        ! Module dependencies
        use biorthogonal_orbital_sets, only: final_orbital_energy_parameters, initial_orbital_energy_parameters, &
            final_orbital_count, initial_orbital_count, final_orbital_principal_n, final_orbital_kappa, &
            initial_orbital_principal_n, initial_orbital_kappa
        use radial_grid_state, only: radial_grid, radial_point_count
        use radial_wavefunction_format, only: write_radial_orbital_set
        use radial_wavefunctions, only: final_orbital_last_grid_index, final_orbital_origin_parameters, &
            initial_orbital_last_grid_index, initial_orbital_origin_parameters, final_radial_large_components, &
            final_radial_small_components, initial_radial_large_components, initial_radial_small_components
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_output_unit
        character(len=128), intent(in) :: name(2)

        open(newunit=wavefunction_output_unit, file=trim(name(1))//".bw", &
            access="stream", form="unformatted", action="write", status="replace")
        write (*, *) "NWII", initial_orbital_count
        call write_radial_orbital_set(wavefunction_output_unit, initial_orbital_count, &
            initial_orbital_principal_n, initial_orbital_kappa, initial_orbital_energy_parameters, &
            initial_orbital_last_grid_index, initial_orbital_origin_parameters, &
            initial_radial_large_components, initial_radial_small_components, radial_grid(:radial_point_count))
        close(wavefunction_output_unit)

        open(newunit=wavefunction_output_unit, file=trim(name(2))//".bw", &
            access="stream", form="unformatted", action="write", status="replace")
        write (*, *) "NWFF", final_orbital_count
        call write_radial_orbital_set(wavefunction_output_unit, final_orbital_count, &
            final_orbital_principal_n, final_orbital_kappa, final_orbital_energy_parameters, &
            final_orbital_last_grid_index, final_orbital_origin_parameters, &
            final_radial_large_components, final_radial_small_components, radial_grid(:radial_point_count))
        close(wavefunction_output_unit)

        return
    end subroutine write_transformed_radial_wavefunctions
end module write_transformed_radial_wavefunctions_interface
