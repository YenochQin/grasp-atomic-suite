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
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:28   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use debug_c
      use def_c
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: ainfcm = 0.52917721067d-08
      real(double), parameter :: alfai  = 137.035999139d00
      real(double), parameter :: ccmps  = 2.99792458d10
      real(double), parameter :: eesu   = 4.803204673d-10
      real(double), parameter :: emeg   = 9.10938356d-28
      real(double), parameter :: emeamu = 5.48579909070d-04
      real(double), parameter :: empamu = 1.007276466879d00
      real(double), parameter :: hbares = 1.054571800d-27
      real(double), parameter :: rinfev = 13.605693009d00
      real(double), parameter :: rinfk  = 109737.31568508d00
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
      empam = empamu
      rbcm = ainfcm
!
!   Calculate constants for  /DEF11/:
!
!      CVAC is the speed of light in the vacuum in atomic units
!      B1 is the conversion factor between the amu and the atomic
!         unit of mass
!
      cvac = alfai
      aumamu = emeamu
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
      aucm = 2.0d00*rinfk
      auev = 2.0d00*rinfev
      ccms = ccmps
      fasi = (emeg/hbares)*(eesu*eesu/hbares)**2
      fbsi = (10.0d00*ainfcm**3/hbares)*fasi
!
!   Calculate conversion factor for Fermis to Bohr radii
!
      fmtoau = 1.0d-13/ainfcm
!
!   Calculate \pi from FORTRAN function
!
      pi = 4.0d00*atan(1.0d00)
!
!   Printouts
!
      if (ldbpg(2)) write (99, 300) ainfcm, alfai, ccmps, eesu, emeg, emeamu, &
         empamu, hbares, rinfev, rinfk
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
