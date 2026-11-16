!***********************************************************************
!                                                                      *
      subroutine ldcsl2(ncore, name)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: IQ, LENGTH, LODCSL, OPENFL.                   *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use orb_c
      use biorb_c

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodcsl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: ncore
      character , intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, ios, ncorer, i
      character :: record*15
!-----------------------------------------------
!
!
!     Common relevant for the radial wave functions of the final state.
!
!
!   The  .csl  file is FORMATTED; it must exist
!
      j = index(name,' ')
      open(unit=21, file=name(1:j-1)//'.c', form='FORMATTED', status='OLD', &
         position='asis')
!
!   Check the first record of the file; if not as expected, try again
!
      read (21, '(1A15)', iostat=ios) record
      if (ios/=0 .or. record(1:15)/='Core subshells:') then
         write (6, *) 'Not a Configuration Symmetry List File;'
         close(21)
         stop
      endif
!
!   Load data from the  .csl  file
!
      call lodcsl (ncorer)
!
!   Close the  .csl  file
!
      close(21)
!
!   Check if the core should be redefined
!
      ncore = ncorer
!      DO 3 I = NCORER+1,NW
!         IFULLI = NKJ(I)+1
!         DO 2 J = 1,NCF
!            IF (IQ (I,J) .NE. IFULLI) GOTO 4
!    2    CONTINUE
!         CALL CONVRT (NP(I),RECORD,LENTH)
!         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
!         PRINT *, ' in all CSFs; including this'
!         PRINT *, ' subshell in the core;'
!         NCORE = NCORE+1
!    3 CONTINUE
!
      nelecff = nelec
      nwff = nw
      ncfff = ncf
      nhff(:nw) = nh(:nw)
      npff(:nw) = np(:nw)
      nakff(:nw) = nak(:nw)
      nklff(:nw) = nkl(:nw)
      nkjff(:nw) = nkj(:nw)
      return
      end subroutine ldcsl2
