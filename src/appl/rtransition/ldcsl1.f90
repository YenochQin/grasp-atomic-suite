!***********************************************************************
!                                                                      *
      subroutine ldcsl1 (ncorer,name)
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, IQ, LENGTH, LODCSL, OPENFL.            *
!                                                                      *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnw
      use memory_man
      use def_c
      use biorb_c
      use orb_c
      use stat_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      use lodcsl_i
      use openfl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: ncorer
      character(len=24) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: i, j, ios, ncore
       character(len=15) :: record
!
!   Check the first record of the file; if not as expected, try again
!
      j = index(name,' ')
      open (unit = 21,file=name(1:j-1)//'.c',form='FORMATTED',          &
     &     status='OLD')

      read (21,'(1A15)',iostat = ios) record
      if ((ios .ne. 0) .or.                                             &
     &    (record(1:15) .ne. 'Core subshells:')) then
         print *, 'Not a Configuration Symmetry List File;'
         close (21)
         stop
      endif
!
!   Load data from the  .csl  file
!
      call lodcsl (ncore)
!
!   Close the  .csl  file
!
      close (21)
!
!   Check if the core should be redefined
!
      ncorer = ncore
!      DO 3 I = NCORE+1,NW
!         IFULLI = NKJ(I)+1
!         DO 2 J = 1,NCF
!            IF (IQ (I,J) .NE. IFULLI) GOTO 4
!    2    CONTINUE
!         CALL CONVRT (NP(I),RECORD,LENTH)
!         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
!         PRINT *, ' in all CSFs; including this'
!         PRINT *, ' subshell in the core;'
!         NCORER = NCORER+1
!    3 CONTINUE
!
!   Copy data to reference arrays
!
    4 nelecr = nelec
      nelecii = nelec
      ncfr = ncf
      nwr = nw
!
      call alloc (iqar, nnnw, ncf, 'IQAR', 'LDCSL1' )
      call alloc (jqsar,nnnw,3,ncf, 'JQSAR', 'LDCSL1')
      call alloc (jcupar,nnnw, ncf, 'JCUPAR', 'LDCSL1')
!
      do 6 i = 1,ncf
         do 5 j = 1,nnnw
            iqar(j,i) = iqa(j,i)
            jqsar(j,1,i) = jqsa(j,1,i)
            jqsar(j,2,i) = jqsa(j,2,i)
            jqsar(j,3,i) = jqsa(j,3,i)
            jcupar(j,i) = jcupa(j,i)
    5    continue
    6 continue
!
      nwii=nw
      ncfii=ncf
      do 7 i = 1,nw
         nhii(i) = nh(i)
         npii(i) = np(i)
         nakii(i) = nak(i)
         nklii(i) = nkl(i)
         nkjii(i) = nkj(i)
         npr(i) = np(i)
         nakr(i) = nak(i)
         nhr(i) = nh(i)
    7 continue
!
!   Deallocate storage so that LODCSL can restart correctly
!
      call dalloc (iqa, 'IQA', 'LDCSL1')
      call dalloc (jqsa, 'JQSA', 'LDCSL1')
      call dalloc (jcupa, 'JCUPA', 'LDCSL1')
!
      close (21)

      return
      end
