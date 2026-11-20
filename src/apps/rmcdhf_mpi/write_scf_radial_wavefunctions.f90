!***********************************************************************
subroutine write_scf_radial_wavefunctions(output_wavefunction_path)
    ! Module dependencies
    !***********************************************************************
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_energy_parameters, orbital_kappa, &
        orbital_principal_n
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    implicit none (type, external)

    ! Arguments and local declarations
    integer :: wavefunction_output_unit
    character , intent(in) :: output_wavefunction_path*(*)
    integer :: orbital_index, last_radial_index, i
    !-----------------------------------------------
    !
    open(newunit=wavefunction_output_unit, file=output_wavefunction_path, status="UNKNOWN", form="UNFORMATTED", position=&
        "asis")
    write (wavefunction_output_unit) "G92RWF"
    do orbital_index = 1, orbital_count
        last_radial_index = orbital_last_grid_index(orbital_index)
        write (wavefunction_output_unit) orbital_principal_n(orbital_index), orbital_kappa(orbital_index), &
            & orbital_energy_parameters(orbital_index), last_radial_index
        write (wavefunction_output_unit) orbital_origin_parameters(orbital_index), (radial_large_components(i,orbital_index), &
            & i=1,last_radial_index), &
            (radial_small_components(i,orbital_index),i=1,last_radial_index)
        write (wavefunction_output_unit) (radial_grid(i),i=1,last_radial_index)               ! This is a waste of resources
    end do
    close(wavefunction_output_unit)

    return
end subroutine write_scf_radial_wavefunctions
