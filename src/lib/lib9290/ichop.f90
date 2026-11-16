!***********************************************************************
!                                                                      *
      integer function ichop (isubsh, icsf)
!                                                                      *
!   ICHOP is -1 if subshell ISUBSH is empty in CSF  ICSF,  +1 if the   *
!   subshell is full, and 0 if it is open.                             *
!                                                                      *
!   Call(s) to: [LIB92]: IQ.                                           *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 30 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:25   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use orb_c,         only: orbital_l, orbital_two_j
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: isubsh
      integer :: icsf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iocc, ifull
!-----------------------------------------------
!
!
! cff  Since ICHOP is always called from within a do-loop over the
!      appropriate range, testing seems redundant
!     IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NW)) THEN
!        IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
      iocc = iq(isubsh,icsf)
      ifull = orbital_two_j(isubsh) + 1
      if (iocc == 0) then
         ichop = -1
      else if (iocc == ifull) then
         ichop = 1
      else
         ichop = 0
      endif
!        ELSE
!           PRINT *, 'ICHOP: Argument ICSF is out of range.'
!           STOP
!        ENDIF
!     ELSE
!        PRINT *, 'ICHOP: Argument ISUBSH is out of range.'
!        STOP
!     ENDIF
!
      return
      end function ichop
