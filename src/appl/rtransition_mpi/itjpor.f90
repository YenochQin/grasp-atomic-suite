!***********************************************************************
!                                                                      *
      integer function itjpor (icsf)
!                                                                      *
!   ITJPOR is the value of 2J+1 for CSF number ICSF.                   *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters, only: max_orbitals
      use csf_coupling_data,        only: jcupar
      use standard_io_units,      only: stderr_unit
      use orbital_csf_basis,         only: ncfr
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: icsf
!-----------------------------------------------
!
      if (icsf>=1 .and. icsf<=ncfr) then
        itjpor = jcupar(max_orbitals,icsf)
         if (itjpor > 127) itjpor = 256 - itjpor
         itjpor = iabs (itjpor)
      else
         write (6, *) "ITJPOR: Argument ICSF is out of range."
         stop
      end if
!
      return
      end function itjpor
