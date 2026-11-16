!***********************************************************************
!                                                                      *
      integer function isparr (icsf)
!                                                                      *
!   ISPARR is the value of P for CSF number ICSF.                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
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
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
      if (icsf>=1 .and. icsf<=ncfr) then
         isparr = jcupar(max_orbitals,icsf)

         if (isparr > 127) isparr = isparr - 256
         isparr = sign(1,isparr)
      else
         write (6, *) "ISPARR: Argument ICSF is out of range."
         stop
      end if
!
      return
      end function isparr
