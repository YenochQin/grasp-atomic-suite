!***********************************************************************
!                                                                      *
      subroutine setcon
!                                                                      *
!   This  subprogram  sets the values of the fundamental and derived   *
!   physical constants, and other useful constants.                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c
      use def_c
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
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
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
!   Physical constants: see description below.
!
!
!   The values of the physical constants used here are taken from
!    (see setcon.f-????)
!   E R Cohen and B N Taylor, The 1986 Adjstment of the Fundamental
!   Physical Constants, Report of the CODATA Task Group on Funda-
!   mental COnstants, CODATA Bulletin 63, Pergamon, Elmsford, NY
!   (1986)
!
!      AINFCM: Bohr radius in cm
!      ALFAI : Inverse of the fine-structure constant
!      CCMPS : Speed of light in cm/s
!      EESU  : Electron charge in esu
!              conversion : EESU = EEC * 2997924580
!      EMEG  : Electron mass in g
!      EMEAMU: Electron mass in amu
!      EMPAMU: Proton mass in amu
!      HBARES: Rationalized Planck's constant in erg s
!      RINFEV: Rydberg in eV
!      RINFK : Rydberg in Kaysers
!
!   In the context of new results, Taylor warns that ... since the
!   output values of a least-squares adjustment are related in a
!   complex way and a change in the measured value of one constant
!   usually leads to corresponding changes in the adjusted values
!   of others, one must be cautious in carrying out calculations
!   using both the [above values] and the results of more recent
!   measurements.
!
!   Calculate constants for /DEF3/:
!
      proton_mass_u = reference_proton_mass_u
      bohr_radius_cm = reference_bohr_radius_cm
!
!   Calculate constants for  /DEF11/:
!
!      CVAC is the speed of light in the vacuum in atomic units
!      B1 is the conversion factor between the amu and the atomic
!         unit of mass
!
      vacuum_speed_of_light_au = inverse_fine_structure_constant
      electron_mass_u = reference_electron_mass_u
!
!   Calculate constants for  /DEF10/:
!
!      AUCM converts au to kaysers assuming an infinitely-heavy
!           nucleus
!      AUEV Converst au to eV assuming an infinitely-heavy
!           nucleus
!      CCMS is the speed of light in the vacuum in centimetres
!            per second
!      FASI converts the Einstein A coefficients from atomic to SI
!            units
!      FBSI converts the Einstein B coefficients from atomic to SI
!            units
!
      hartree_to_wavenumber_cm = 2.0_real64*rydberg_wavenumber_cm
      hartree_to_ev = 2.0_real64*rydberg_energy_ev
      speed_of_light_cm_per_s = reference_speed_of_light_cm_per_s
      einstein_a_au_to_si = (electron_mass_g/reduced_planck_constant_erg_s)*(electron_charge_esu*electron_charge_esu/ &
         & reduced_planck_constant_erg_s)**2
      einstein_b_au_to_si = (10.0_real64*reference_bohr_radius_cm**3/reduced_planck_constant_erg_s)*einstein_a_au_to_si
!
!   Calculate conversion factor for Fermis to Bohr radii
!
      femtometer_to_bohr = 1.0e-13_real64/reference_bohr_radius_cm
!
!   Calculate \pi from FORTRAN function
!
      pi = 4.0_real64*atan(1.0_real64)
!
!   Printouts
!
      if (ldbpg(2)) write (99, 300) reference_bohr_radius_cm, inverse_fine_structure_constant, &
         & reference_speed_of_light_cm_per_s, electron_charge_esu, electron_mass_g, reference_electron_mass_u, &
         reference_proton_mass_u, reduced_planck_constant_erg_s, rydberg_energy_ev, rydberg_wavenumber_cm
!
      return
!
  300 format(/,'From SUBROUTINE SETCON:'/,' AINFCM (Bohr radius in cm): ',0p,1d&
         16.9,','/,' ALFAI (Inverse of the fine-structure constant): ',3p,1d&
         16.9,','/,' CCMPS (Speed of light in cm/s): ',1p,1d15.8,','/,&
         ' EESU (Electron charge in esu): ',1p,1d15.8,','/,&
         ' EMEG (Electron mass in g): ',1p,1d14.7,','/,&
         ' EMEAMU (Electron mass in u): ',1p,1d16.9,','/,&
         ' EMPAMU (Proton mass in u): ',1p,1d16.9,','/,&
         ' HBARES (Rationalized Planck constant in erg s): ',1p,1d15.8,','/,&
         ' RINFEV (Rydberg in eV): ',2p,1d15.8,','/,&
         ' RINFK (Rydberg in Kaysers): ',6p,1d17.10,'.')
      return
!
      end subroutine setcon
