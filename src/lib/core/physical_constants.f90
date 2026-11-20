! Initialize physical constants and conversions between atomic and conventional units.
module physical_constants
    implicit none (type, external)
    private
    public :: setcon
contains
      subroutine setcon
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: bohr_radius_cm, einstein_a_au_to_si, einstein_b_au_to_si, electron_mass_u, &
          femtometer_to_bohr, hartree_to_ev, hartree_to_wavenumber_cm, pi, proton_mass_u, speed_of_light_cm_per_s, &
          vacuum_speed_of_light_au
      implicit none (type, external)
      real(real64), parameter :: reference_bohr_radius_cm = 0.52917721067e-8_real64
      real(real64), parameter :: inverse_fine_structure_constant  = 137.035999139_real64
      real(real64), parameter :: reference_speed_of_light_cm_per_s  = 2.99792458e10_real64
      real(real64), parameter :: electron_charge_esu   = 4.803204673e-10_real64
      real(real64), parameter :: electron_mass_g   = 9.10938356e-28_real64
      real(real64), parameter :: reference_electron_mass_u = 5.48579909070e-4_real64
      real(real64), parameter :: reference_proton_mass_u = 1.007276466879_real64
      real(real64), parameter :: reduced_planck_constant_erg_s = 1.054571800e-27_real64
      real(real64), parameter :: rydberg_energy_ev = 13.605693009_real64
      real(real64), parameter :: rydberg_wavenumber_cm  = 109737.31568508_real64
      proton_mass_u = reference_proton_mass_u
      bohr_radius_cm = reference_bohr_radius_cm
      vacuum_speed_of_light_au = inverse_fine_structure_constant
      electron_mass_u = reference_electron_mass_u
      hartree_to_wavenumber_cm = 2.0_real64*rydberg_wavenumber_cm
      hartree_to_ev = 2.0_real64*rydberg_energy_ev
      speed_of_light_cm_per_s = reference_speed_of_light_cm_per_s
      einstein_a_au_to_si = (electron_mass_g/reduced_planck_constant_erg_s)*(electron_charge_esu*electron_charge_esu/ &
         & reduced_planck_constant_erg_s)**2
      einstein_b_au_to_si = (10.0_real64*reference_bohr_radius_cm**3/reduced_planck_constant_erg_s)*einstein_a_au_to_si
      femtometer_to_bohr = 1.0e-13_real64/reference_bohr_radius_cm
      pi = 4.0_real64*atan(1.0_real64)
      if (ldbpg(2)) write (99, '(/,"From SUBROUTINE SETCON:"/," AINFCM (Bohr radius in cm): ",0p,1d16.9,","/,&
         &" ALFAI (Inverse of the fine-structure constant): ",3p,1d16.9,","/," CCMPS (Speed of light in cm/s): ",&
         &1p,1d15.8,","/," EESU (Electron charge in esu): ",1p,1d15.8,","/," EMEG (Electron mass in g): ",1p,&
         &1d14.7,","/," EMEAMU (Electron mass in u): ",1p,1d16.9,","/," EMPAMU (Proton mass in u): ",1p,1d16.9,","/&
         &," HBARES (Rationalized Planck constant in erg s): ",1p,1d15.8,","/," RINFEV (Rydberg in eV): ",2p,&
         &1d15.8,","/," RINFK (Rydberg in Kaysers): ",6p,1d17.10,".")') &
         reference_bohr_radius_cm, inverse_fine_structure_constant, &
         & reference_speed_of_light_cm_per_s, electron_charge_esu, electron_mass_g, reference_electron_mass_u, &
         reference_proton_mass_u, reduced_planck_constant_erg_s, rydberg_energy_ev, rydberg_wavenumber_cm
      return
      return
      end subroutine setcon
end module physical_constants
