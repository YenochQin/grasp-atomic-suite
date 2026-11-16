!***********************************************************************
!                                                                      *
      subroutine strsum(name, inpci, ilbl)
!                                                                      *
!   Generates the first part of  oscl92.sum  (on stream 24).           *
!                                                                      *
!   Call(s) to: [LIB92]: CALEN, CONVRT, WGHTD5.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use biorb_c
      use decide_c
      use def_c
      use foparm_c
      use grid_c
      use eigv_c
      use npar_c
      use orb_c
      use prnt_c
      use syma_c
      use wave_c
      use jj2lsjbio_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: inpci
      integer :: ilbl
      character, intent(in) :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      character :: record*15, ctime*8, cdate*8
!-----------------------------------------------


      i = index(name(1),' ')
      j = index(name(2),' ')
      if(ilbl == 0) then
         if (inpci == 0) then
            open(unit=24,file=name(1)(1:i-1)//'.'//name(2)(1:j-1)//'.ct',form=&
            'FORMATTED', status='UNKNOWN',position='asis')
         else
            open(unit=24,file=name(1)(1:i-1)//'.'//name(2)(1:j-1)//'.t',form=&
            'FORMATTED',status='UNKNOWN',position='asis')
         endif
      else if(ilbl == 1) then
         if(iopen_status1.eq.0 .and. iopen_status2 .eq.0) then
            if (inpci == 0) then
               open(unit=32,                                           &
                   file=name(1)(1:i-1)//'.'//name(2)(1:j-1)//'.ct.lsj',&
                   form='FORMATTED',status='UNKNOWN')
            else
               open(unit=32,                                           &
                   file=name(1)(1:i-1)//'.'//name(2)(1:j-1)//'.t.lsj', &
                   form='FORMATTED',status='UNKNOWN')
            endif
         end if
      end if
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!ww      CALL CALEN (CTIME,CDATE)
!ww      WRITE (24,*) 'OSCL92 run at ',CTIME,' on ',CDATE,'.'
!
!   Write out the basic dimensions of the initial state electron cloud
!
!ww      WRITE (24,*)
!ww      CALL CONVRT (NELECII,RECORD,LENTH)
!ww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
!ww     :           //' electrons in the initial state cloud'
!ww      CALL CONVRT (NCFII,RECORD,LENTH)
!ww      WRITE (24,*) ' in '//RECORD(1:LENTH)
!ww     :           //' relativistic CSFs'
!ww      CALL CONVRT (NWII,RECORD,LENTH)
!ww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
!ww     :           //' relativistic subshells.'
!
!   Write out the basic dimensions of the final state electron cloud
!
!ww      WRITE (24,*)
!ww      CALL CONVRT (NELECFF,RECORD,LENTH)
!ww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
!ww     :           //' electrons in the final state cloud'
!ww      CALL CONVRT (NCFFF,RECORD,LENTH)
!ww      WRITE (24,*) ' in '//RECORD(1:LENTH)
!ww     :           //' relativistic CSFs'
!ww      CALL CONVRT (NWFF,RECORD,LENTH)
!ww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
!ww     :           //' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
!ww      IF (LFORDR) THEN
!ww         WRITE (24,*)
!ww         CALL CONVRT (ICCUT,RECORD,LENTH)
!ww         WRITE (24,*) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'
!ww     :              //' the zero-order space;'
!ww      ENDIF
!
!   Write out the nuclear parameters
!
!ww      WRITE (24,*)
!ww      WRITE (24,300) Z
!ww      IF (NPARM .EQ. 2) THEN
!ww         WRITE (24,*) 'Fermi nucleus:'
!ww         WRITE (24,301) PARM(1),PARM(2)
!ww      ELSE
!ww         WRITE (24,*) ' point nucleus.'
!ww      ENDIF
!
!   Write out the physical effects specifications
!
!ww      WRITE (24,*)
!ww      WRITE (24,305) C
!
!   Write out the parameters of the radial grid
!
!ww      WRITE (24,*)
!ww      IF (HP .EQ. 0.0D 00) THEN
!ww         WRITE (24,306) RNT,H,N
!ww      ELSE
!ww         WRITE (24,307) RNT,H,HP,N
!ww      ENDIF
!ww      WRITE (24,308) R(1),R(2),R(N)
!
!   Write out the orbital properties
!
!ww      WRITE (24,*)
!ww      WRITE (24,*) 'Initial state subshell radial wavefunction summary:'
!ww      WRITE (24,*)
!ww      WRITE (24,309)
!ww      WRITE (24,*)
!ww      DO 1 I = 1,NWII
!ww         WRITE (24,310) NPII(I),NHII(I),EII(I),PZII(I),
!ww     :                  GAMAII(I),PFII(2,I),QFII(2,I),MFII(I)
!ww    1 CONTINUE
!
!ww      WRITE (24,*)
!ww      WRITE (24,*) 'Final state subshell radial wavefunction summary:'
!ww      WRITE (24,*)
!ww      WRITE (24,309)
!ww      WRITE (24,*)
!ww      DO 2 I = 1,NWFF
!ww         WRITE (24,310) NPFF(I),NHFF(I),EFF(I),PZFF(I),
!ww     :                  GAMAFF(I),PFFF(2,I),QFFF(2,I),MFFF(I)
!ww    2 CONTINUE
!
!   Write the list of eigenpair indices for the initial state
!
!     WRITE (24,*)
!     CALL ENGOUT1 (EAVII,EVALII,IATJPOII,IASPARII,IVECII,NVECII,3,1)
!
!   Write the list of eigenpair indices for the final state
!
!     WRITE (24,*)
!     CALL ENGOUT1 (EAVFF,EVALFF,IATJPOFF,IASPARFF,IVECFF,NVECFF,3,2)
!
      return
!
  300 format ('The atomic number is ',1f14.10,';')
  301 format ('  c =',1p,1d19.12,' Bohr radii,'            &
             /'  a =',   1d19.12,' Bohr radii;')
  305 format ('Speed of light = ',1pd19.12,' atomic units.')
  306 format ( 'Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',&
               ' I = 1, ..., N;'                           &
             //' RNT  = ',1p,d19.12,' Bohr radii;'         &
              /' H    = ',   d19.12,' Bohr radii;'         &
              /' N    = ',1i4,';')
  307 format ( 'Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',&
               ' I = 1, ..., N;'                           &
             //' RNT  = ',1p,d19.12,' Bohr radii;'         &
              /' H    = ',   d19.12,' Bohr radii;'         &
              /' HP   = ',   d19.12,' Bohr radii;'         &
              /' N    = ',1i4,';')
  308 format ( ' R(1) = ',1p,1d19.12,' Bohr radii;'        &
              /' R(2) = ',   1d19.12,' Bohr radii;'        &
              /' R(N) = ',   1d19.12,' Bohr radii.')
  309 format (' Subshell',11x,'e',20x,'p0',18x,            &
              'gamma',19x,'P(2)',18x,'Q(2)',10x,'MTP')
  310 format (3x,1i2,1a2,1x,1p,5(3x,1d19.12),3x,1i3)

      return
!
      end subroutine strsum
