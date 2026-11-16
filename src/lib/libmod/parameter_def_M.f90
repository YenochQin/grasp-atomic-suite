!====================================================================
      module parameter_def
!    THis module defines some global parameters for the
!    application
!...Translated by Pacific-Sierra Research 77to90 4.3E 11:01:42 1/15/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!====================================================================
!    KEYORB: Is a packing parameter defined by a packing algorithm.
!            This is the largest integer for which KEYORB**4 < 2^31
!            (the largest positive  INTEGER*4)
!    NNNP:   Number  of points in the radial grid
!    NNN1:   = NNNP+10
!    NNNW    Maximum number or orbitals (previously 120). An n=10
!            calculation has 100 orbitals.  This parameter is used
!            to assign array dimensions such as NNNW*NCFG where
!            NCFG is the total number of CSFs.
!    NNNWM1: = NNNW-1
!    NNNWM2: = NNNW-2
!------------------------------------------------------------------

!    Compatibility interface: all limits are defined in suite_parameters.
!    Keep existing USE parameter_def callers without duplicate definitions.
      use suite_parameters, only: orbital_index_radix
      use suite_parameters, only: max_radial_points
      use suite_parameters, only: radial_workspace_capacity
      use suite_parameters, only: max_orbitals
      use suite_parameters, only: max_orbitals_minus_one
      use suite_parameters, only: max_orbitals_minus_two
! Legacy capacity names remain aliases of the same named constants.
      use suite_parameters, only: keyorb => orbital_index_radix
      use suite_parameters, only: nnnp => max_radial_points
      use suite_parameters, only: nnn1 => radial_workspace_capacity
      use suite_parameters, only: nnnw => max_orbitals
      use suite_parameters, only: nnnwm1 => max_orbitals_minus_one
      use suite_parameters, only: nnnwm2 => max_orbitals_minus_two
      implicit none
      end module parameter_def
!======================================================================
