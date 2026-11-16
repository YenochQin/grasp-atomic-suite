!***********************************************************************
!                                                                      *
      subroutine lodiso
!                                                                      *
!   Loads the data from the  .iso  file.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 29 Sep 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use def_c
      use npar_c
      use nsmdat_c,        only: sqn, dmomnm, qmomb
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: a, aparm, cparm, emnamu
!
!   Read and echo pertinent information from  .iso  file
!
!   Atomic number
!
      read (22, *) nuclear_charge
!
!   Nuclear geometry
!
      read (22, *)
      read (22, *) a
      read (22, *)
      read (22, *) aparm
      read (22, *)
      read (22, *) cparm
!
      if (a /= 0.d0) then
         nuclear_model_parameter_count = 2
         nuclear_model_parameters(1) = cparm*femtometer_to_bohr
         nuclear_model_parameters(2) = aparm*femtometer_to_bohr
      else
         nuclear_model_parameter_count = 0
      endif
!
!   Nuclear mass
!
      read (22, *)
      read (22, *) emnamu

!
      if (emnamu /= 0.d0) then
         emn = emnamu/electron_mass_u
      else
         emn = 0.d0
      endif


!
!   Nuclear spin and moments
!
      read (22, *)
      read (22, *) sqn
      read (22, *)
      read (22, *) dmomnm
      read (22, *)
      read (22, *) qmomb
!
!   Grid parameters from isodata
!   Jon Grumer (Lund, 2013)
!
!     READ (22,*)
!     READ (22,*) RNT
!     READ (22,*)
!     READ (22,*) H
!     READ (22,*)
!     READ (22,*) HP
!     READ (22,*)
!     READ (22,*) N

      return
      end subroutine lodiso
