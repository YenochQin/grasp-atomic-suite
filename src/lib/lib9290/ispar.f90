!***********************************************************************
!                                                                      *
      integer function ispar (icsf)
!                                                                      *
!   ISPAR is the value of P for CSF number ICSF.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:41   2/14/04
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
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
      if (icsf>=1 .and. icsf<=ncf) then
         ispar = jcupa(nnnw,icsf)
         if (ispar > 127) ispar = ispar - 256
         ispar = sign(1,ispar)
      else
         write (istde, *) 'ISPAR: Argument ICSF is out of range.'
         write (istde, *) 'ICSF =',icsf,' NCF =',ncf
         stop
      endif
!
      return
      end function ispar
