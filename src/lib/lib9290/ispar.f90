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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use stat_c,        only: csf_coupling_quantum_codes
      use iounit_c,      only: stderr_unit
      use orb_c,         only: csf_count
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: icsf
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
      if (icsf>=1 .and. icsf<=csf_count) then
         ispar = csf_coupling_quantum_codes(max_orbitals,icsf)
         if (ispar > 127) ispar = ispar - 256
         ispar = sign(1,ispar)
      else
         write (stderr_unit, *) 'ISPAR: Argument ICSF is out of range.'
         write (stderr_unit, *) 'ICSF =',icsf,' NCF =',csf_count
         stop
      endif
!
      return
      end function ispar
