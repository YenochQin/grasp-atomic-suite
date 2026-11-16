module jcup_i
    implicit none (type, external)
    private
    public :: jcup
contains
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use iounit_c,      only: stderr_unit
      use orb_c,         only: orbital_count, csf_count
      use stat_c,        only: csf_coupling_quantum_codes
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: loc
      integer, intent(in) :: icsf
!-----------------------------------------------
!
      if (loc>=1 .and. loc<=max_orbitals-1) then
         if (icsf>=1 .and. icsf<=csf_count) then
            jcup = csf_coupling_quantum_codes(loc,icsf)
         else
            write (stderr_unit, *) 'JCUP: Argument ICSF is out of range.'
            stop
         endif
      else
         write (stderr_unit, *) 'JCUP: Argument LOC is out of range.'
         stop
      endif
!
      return
      end function jcup
end module jcup_i
