module write_initial_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_initial_radial_wavefunctions
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_initial_radial_wavefunctions(file_path)
        ! Module dependencies
        !                                                                      *
        !   Open, write a header and all subshell radial wavefunctions, and    *
        !   close the  .rwf  file.                                             *
        !                                                                      *
        !   Call(s) to: [LIB92]: DALLOC, OPENFL.                               *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        !
        use array_memory_management, only: dalloc
        use radial_grid_state, only: radial_grid
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_kappa, orbital_principal_n
        use radial_wavefunctions, only: orbital_origin_parameters, radial_large_components, radial_small_components, &
            orbital_last_grid_index
        use standard_io_units, only: stderr_unit
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in) :: file_path

        ! Local declarations
        integer :: wavefunction_output_unit, open_status, orbital_index, last_radial_index, radial_index
        !-----------------------------------------------
        !
        !


        call open_calculation_file (wavefunction_output_unit, file_path, "UNFORMATTED", "NEW", open_status)
        if (open_status == 1) then
            write (stderr_unit, *) 'Error when opening "', file_path(1:len_trim(file_path)), &
                '"'

            stop
        end if
        !
        !   Write the file header
        !
        write (wavefunction_output_unit) "G92RWF"
        !
        !   Write out the radial wavefunctions
        !
        do orbital_index = 1, orbital_count
            last_radial_index = orbital_last_grid_index(orbital_index)
            write (wavefunction_output_unit) orbital_principal_n(orbital_index), orbital_kappa(orbital_index), &
                & orbital_energy_parameters(orbital_index), last_radial_index
            write (wavefunction_output_unit) orbital_origin_parameters(orbital_index), &
                & (radial_large_components(radial_index,orbital_index), &
                & radial_index=1,last_radial_index), &
                & (radial_small_components(radial_index,orbital_index),radial_index=1,last_radial_index)
            write (wavefunction_output_unit) (radial_grid(radial_index),radial_index=1,last_radial_index)
        end do

        close(wavefunction_output_unit)
        !
        !   Deallocate the storage for the radial wavefunctions
        !
        call dalloc (radial_large_components, "radial_large_components", "write_initial_radial_wavefunctions")
        call dalloc (radial_small_components, "radial_small_components", "write_initial_radial_wavefunctions")
        !
        return
    end subroutine write_initial_radial_wavefunctions
end module write_initial_radial_wavefunctions_interface
