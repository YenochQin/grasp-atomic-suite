!***********************************************************************
!                                                                      *
      integer function jcup (loc, icsf)
!                                                                      *
!   JCUP is the 2J+1 value of the LOCth nontrivial intermediate ang-   *
!   ular momentum in CSF  ICSF.                                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  08:16:51   2/21/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: nnnw
      use iounit_c,      only: istde
      use orb_c,         only: nw, ncf
      use stat_c,        only: jcupa
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: loc
      integer             :: icsf
!-----------------------------------------------
!
      if (loc>=1 .and. loc<=nnnw-1) then
         if (icsf>=1 .and. icsf<=ncf) then
            jcup = jcupa(loc,icsf)
         else
            write (istde, *) 'JCUP: Argument ICSF is out of range.'
            stop
         endif
      else
         write (istde, *) 'JCUP: Argument LOC is out of range.'
         stop
      endif
!
      return
      end function jcup
