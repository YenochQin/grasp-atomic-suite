!***********************************************************************
!                                                                      *
      integer function itjpo (icsf)
!                                                                      *
!   ITJPO is the value of 2J+1 for CSF number ICSF.                    *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:45   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: nnnw
      use stat_c,        only: jcupa
      use iounit_c,      only: istde
      use orb_c,         only: ncf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer       :: icsf
!-----------------------------------------------
      if (icsf>=1 .and. icsf<=ncf) then
        itjpo = jcupa(nnnw,icsf)
        if (itjpo > 127) itjpo = 256 - itjpo
        itjpo = iabs (itjpo)
      else
         write (istde, *) 'ITJPO: Argument ICSF is out of range.'
         stop
      endif
!
      return
      end function itjpo
