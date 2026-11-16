!***********************************************************************
!                                                                      *
      integer function jcupr (loc, icsf)
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
      use parameter_def, only: max_orbitals
      use orb_c,         only: nwr, ncfr
      use stat_c,        only: jcupar
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: loc
      integer, intent(in) :: icsf
!-----------------------------------------------
!
      if (loc>=1 .and. loc<=nwr-1) then
         if (icsf>=1 .and. icsf<=ncfr) then
            jcupr = jcupar(loc,icsf)
         else
            write (6, *) 'JCUPR: Argument ICSF is out of range.'
            stop
         endif
      else
         write (6, *) 'JCUPR: Argument LOC is out of range.'
         stop
      endif
!
      return
      end function jcupr
