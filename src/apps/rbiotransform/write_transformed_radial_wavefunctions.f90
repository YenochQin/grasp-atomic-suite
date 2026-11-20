module write_transformed_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_transformed_radial_wavefunctions
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_transformed_radial_wavefunctions(name)
        ! Module dependencies
        !                                                                      *
        !  This subroutine outputs the transformed radial orbitals             *
        !                                                                      *
        !  Written by Per Jonsson                                              *
        !                                                                      *
        !***********************************************************************
        use biorthogonal_orbital_sets, only: final_orbital_energy_parameters, initial_orbital_energy_parameters, &
            final_orbital_count, initial_orbital_count, final_orbital_principal_n, final_orbital_kappa, &
            initial_orbital_principal_n, initial_orbital_kappa
        use radial_grid_state, only: radial_grid
        use radial_wavefunctions, only: final_orbital_last_grid_index, final_orbital_origin_parameters, &
            initial_orbital_last_grid_index, initial_orbital_origin_parameters, final_radial_large_components, &
            final_radial_small_components, initial_radial_large_components, initial_radial_small_components
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_output_unit
        character, intent(in) :: name(2)*24
        integer :: j, k, i
        !-----------------------------------------------
        !
        j = index(name(1)," ")
        open(newunit=wavefunction_output_unit, file=name(1)(1:j-1)//".bw", form="UNFORMATTED", status=&
            "UNKNOWN", position="asis")

        write (wavefunction_output_unit) "G92RWF"
        write (*, *) "NWII", initial_orbital_count
        do k = 1, initial_orbital_count
            write (wavefunction_output_unit) initial_orbital_principal_n(k), initial_orbital_kappa(k), &
                initial_orbital_energy_parameters(k), initial_orbital_last_grid_index(k)
            write (wavefunction_output_unit) initial_orbital_origin_parameters(k), (initial_radial_large_components(i,k),i=1, &
                initial_orbital_last_grid_index(k)), (initial_radial_small_components(i,k),i=1,initial_orbital_last_grid_index(k))
            write (wavefunction_output_unit) (radial_grid(i),i=1,initial_orbital_last_grid_index(k))
        end do

        close(wavefunction_output_unit)

        j = index(name(2)," ")
        open(newunit=wavefunction_output_unit, file=name(2)(1:j-1)//".bw", form="UNFORMATTED", status=&
            "UNKNOWN", position="asis")

        write (wavefunction_output_unit) "G92RWF"
        write (*, *) "NWFF", final_orbital_count
        do k = 1, final_orbital_count
            write (wavefunction_output_unit) final_orbital_principal_n(k), final_orbital_kappa(k), &
                final_orbital_energy_parameters(k), final_orbital_last_grid_index(k)
            write (wavefunction_output_unit) final_orbital_origin_parameters(k),(final_radial_large_components(i,k),i=1, &
                final_orbital_last_grid_index(k)),(final_radial_small_components(i,k),i=1,final_orbital_last_grid_index(k))
            write (wavefunction_output_unit) (radial_grid(i),i=1,final_orbital_last_grid_index(k))
        end do

        close(wavefunction_output_unit)

        return
    end subroutine write_transformed_radial_wavefunctions
end module write_transformed_radial_wavefunctions_interface
