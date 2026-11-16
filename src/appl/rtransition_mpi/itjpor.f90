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
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use stat_c,        only: jcupar
      use iounit_c,      only: stderr_unit
      use orb_c,         only: ncfr
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
         write (6, *) 'ITJPOR: Argument ICSF is out of range.'
         stop
      endif
!
      return
      end function itjpor
