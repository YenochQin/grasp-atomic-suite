! Load the nuclear model, mass, spin, and moments from the isotope records on unit 22.
module isotope_data_records
    implicit none (type, external)

    ! Module interface
    private
    public :: lodiso
contains
    subroutine lodiso(isotope_input_unit)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: &
            electron_mass_u, nuclear_mass_in_electron_masses => emn, femtometer_to_bohr, nuclear_charge
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use isotope_properties, only: nuclear_spin => sqn, nuclear_magnetic_dipole_moment => dmomnm, &
            nuclear_electric_quadrupole_moment => qmomb
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: isotope_input_unit

        ! Local declarations
        real(real64) :: nuclear_mass_number, fermi_diffuseness_fm, fermi_half_density_radius_fm, nuclear_mass_u
        read (isotope_input_unit, *) nuclear_charge
        read (isotope_input_unit, *)
        read (isotope_input_unit, *) nuclear_mass_number
        read (isotope_input_unit, *)
        read (isotope_input_unit, *) fermi_diffuseness_fm
        read (isotope_input_unit, *)
        read (isotope_input_unit, *) fermi_half_density_radius_fm
        if (nuclear_mass_number /= 0.0_real64) then
            nuclear_model_parameter_count = 2
            nuclear_model_parameters(1) = fermi_half_density_radius_fm*femtometer_to_bohr
            nuclear_model_parameters(2) = fermi_diffuseness_fm*femtometer_to_bohr
        else
            nuclear_model_parameter_count = 0
        end if
        read (isotope_input_unit, *)
        read (isotope_input_unit, *) nuclear_mass_u

        if (nuclear_mass_u /= 0.0_real64) then
            nuclear_mass_in_electron_masses = nuclear_mass_u/electron_mass_u
        else
            nuclear_mass_in_electron_masses = 0.0_real64
        end if

        read (isotope_input_unit, *)
        read (isotope_input_unit, *) nuclear_spin
        read (isotope_input_unit, *)
        read (isotope_input_unit, *) nuclear_magnetic_dipole_moment
        read (isotope_input_unit, *)
        read (isotope_input_unit, *) nuclear_electric_quadrupole_moment

        return
    end subroutine lodiso
end module isotope_data_records
