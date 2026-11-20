!***********************************************************************
!                                                                      *
      integer function jcupr (loc, csf_index)
      use, intrinsic :: iso_fortran_env, only: output_unit
!                                                                      *
!   JCUPR is the 2J+1 value of the LOCth nontrivial intermediate ang-  *
!   ular momentum in CSF  ICSF.                                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters, only: max_orbitals
      use orbital_csf_basis,         only: nwr, ncfr
      use csf_coupling_data,        only: jcupar
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: loc
      integer, intent(in) :: csf_index
!-----------------------------------------------
!
      if (loc>=1 .and. loc<=nwr-1) then
         if (csf_index>=1 .and. csf_index<=ncfr) then
            jcupr = jcupar(loc,csf_index)
         else
            write (output_unit, *) "JCUPR: Argument ICSF is out of range."
            stop
         end if
      else
         write (output_unit, *) "JCUPR: Argument LOC is out of range."
         stop
      end if
!
      return
      end function jcupr
