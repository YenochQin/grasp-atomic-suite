! Convert nuclear moments for the MHz interaction constants and the Hartree matrix elements.
module hyperfine_nuclear_factors
    implicit none (type, external)
    private
    public :: calculate_nuclear_factors
contains
    subroutine calculate_nuclear_factors(dipole_mhz, quadrupole_mhz, reduced_dipole, reduced_quadrupole)
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: bohr_radius_cm, electron_mass_u, hartree_to_wavenumber_cm, &
            proton_mass_u, speed_of_light_cm_per_s, vacuum_speed_of_light_au
        use isotope_properties, only: nuclear_spin, nuclear_magnetic_dipole_moment, nuclear_electric_quadrupole_moment
        implicit none (type, external)
        real(real64), intent(out) :: dipole_mhz, quadrupole_mhz, reduced_dipole, reduced_quadrupole
        real(real64) :: hartree_to_mhz, barn_to_bohr_squared, nuclear_magneton_au

        dipole_mhz = 0.0_real64
        quadrupole_mhz = 0.0_real64
        reduced_dipole = 0.0_real64
        reduced_quadrupole = 0.0_real64
        ! A spin-zero nucleus has no magnetic dipole or electric quadrupole tensor.
        if (nuclear_spin <= 0.0_real64) then
            if (nuclear_spin < 0.0_real64) error stop "HYPERFINE requires a nonnegative nuclear spin"
            return
        end if

        hartree_to_mhz = hartree_to_wavenumber_cm*speed_of_light_cm_per_s*1.0e-6_real64
        barn_to_bohr_squared = 1.0e-24_real64/bohr_radius_cm**2
        nuclear_magneton_au = electron_mass_u/(2.0_real64*vacuum_speed_of_light_au*proton_mass_u)
        dipole_mhz = hartree_to_mhz*nuclear_magneton_au*nuclear_magnetic_dipole_moment/nuclear_spin
        quadrupole_mhz = hartree_to_mhz*2.0_real64*nuclear_electric_quadrupole_moment*barn_to_bohr_squared
        reduced_dipole = sqrt((nuclear_spin + 1.0_real64)/nuclear_spin)*nuclear_magnetic_dipole_moment* &
            0.5_real64*electron_mass_u/(proton_mass_u*vacuum_speed_of_light_au)
        if (nuclear_spin >= 1.0_real64) then
            reduced_quadrupole = sqrt((nuclear_spin + 1.0_real64)*(2.0_real64*nuclear_spin + 3.0_real64)/ &
                ((2.0_real64*nuclear_spin - 1.0_real64)*nuclear_spin))
            reduced_quadrupole = reduced_quadrupole*nuclear_electric_quadrupole_moment*0.5_real64*barn_to_bohr_squared
        end if
    end subroutine calculate_nuclear_factors
end module hyperfine_nuclear_factors
