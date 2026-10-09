!***********************************************************************
!
      SUBROUTINE GETMIXBLOCK(NAME, NCI)
!
!    Reads mixing coefficient file from block-structured format
!
! Note:
!     eav is not compatible with the non-block version if some blocks
!     were not diagonalized
!
!     This is a modified version of cvtmix.f
!     written by Per Jonsson, September 2003
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:32:57   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/01/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE density_parallel, ONLY: density_fail
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
      USE, INTRINSIC :: iso_fortran_env, ONLY: int64
      USE def_C
      USE EIGV_C
      USE orb_C
      USE prnt_C
      USE syma_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)           :: NCI
      CHARACTER(LEN=24), INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IERR, IOS, NCFTOT, NVECTOT, NVECSIZ, NBLOCK, I, NVECPAT, &
         NCFPAT, NVECSIZPAT, NEAVSUM, JB, NB, NCFBLK, NEVBLK, IATJP, IASPA, J
      INTEGER :: NELECMIX, NWMIX
      REAL(DOUBLE) :: EAVSUM
      CHARACTER :: FILNAM*256, FORM*11, G92MIX*6, STATUS*3
!-----------------------------------------------
!
!   The  .mix  file is UNFORMATTED; it must exist
!
      K = INDEX(NAME,' ')
      IF (NCI == 0) THEN
         FILNAM = NAME(1:K-1)//'.cm'
      ELSE
         FILNAM = NAME(1:K-1)//'.m'
      ENDIF
      FORM = 'UNFORMATTED'
      STATUS = 'OLD'
!
      CALL OPENFL (25, FILNAM, FORM, STATUS, IERR)
      IF (IERR == 1) THEN
         WRITE (ISTDE, *) 'Error when opening', FILNAM
         CALL density_fail('Cannot read mixing coefficients')
      ENDIF
!
!   Check the header of the file; if not as expected, try again
!
      READ (25, IOSTAT=IOS) G92MIX
      IF (IOS/=0 .OR. G92MIX/='G92MIX') THEN
         WRITE (ISTDE, *) 'Not a GRASP92 MIXing Coefficients File;'
         CLOSE(25)
         CALL density_fail('Cannot read mixing coefficients')
      ENDIF

      READ (25, IOSTAT=IOS) NELECMIX, NCFTOT, NWMIX, NVECTOT, NVECSIZ, NBLOCK
      IF (IOS /= 0) CALL density_fail('Truncated mixing header')
      IF (NELECMIX /= NELEC .OR. NCFTOT /= NCF .OR. NWMIX /= NW) &
         CALL density_fail('CSF and mixing-file dimensions do not match')
      IF (NCFTOT < 1 .OR. NVECTOT < 1 .OR. NBLOCK < 1) &
         CALL density_fail('Empty or invalid mixing file')
      IF (INT(NCFTOT,int64)*NVECTOT > HUGE(NCF)) &
         CALL density_fail('Mixing vector size exceeds integer index capacity')
      WRITE (*, *) '   nelec  = ', NELEC
      WRITE (*, *) '   ncftot = ', NCFTOT
      WRITE (*, *) '   nw     = ', NW
      WRITE (*, *) '   nblock = ', NBLOCK
      WRITE (*, *)

!***********************************************************************
! Allocate memory for old format data
!***********************************************************************

      CALL ALLOC (EVAL, NVECTOT, 'EVAL', 'GETMIXBLOCK')
      CALL ALLOC (EVEC, NCFTOT*NVECTOT, 'EVEC', 'GETMIXBLOCK')
      CALL ALLOC (IVEC, NVECTOT, 'IVEC', 'GETMIXBLOCK')
      CALL ALLOC (IATJPO, NVECTOT, 'IATJPO', 'GETMIXBLOCK')
      CALL ALLOC (IASPAR, NVECTOT, 'IASPAR', 'GETMIXBLOCK')

!***********************************************************************
! Initialize mixing coefficients to zero; others are fine
!***********************************************************************
      EVEC(:NVECTOT*NCFTOT) = 0.D0

!***********************************************************************
! Initialize counters and sum registers
!
!    nvecpat:    total number of eigenstates of the previous blocks
!    ncfpat:     total number of CSF of the previous blocks
!    nvecsizpat: vector size of the previous blocks
!    eavsum:     sum of diagonal elements of the previous blocks where
!                at least one eigenstate is calculated
!    neavsum:    total number CSF contributing to eavsum
!***********************************************************************

      NVECPAT = 0
      NCFPAT = 0
      NVECSIZPAT = 0
      NEAVSUM = 0
      EAVSUM = 0.D0

      WRITE (*, *) '  block     ncf     nev    2j+1  parity'
      DO JB = 1, NBLOCK

         READ (25, IOSTAT=IOS) NB, NCFBLK, NEVBLK, IATJP, IASPA
         IF (IOS /= 0) CALL density_fail('Truncated mixing block header')
         IF (NCFBLK < 0 .OR. NEVBLK < 0 .OR. NEVBLK > NCFBLK .OR. &
             NCFBLK > NCFTOT-NCFPAT .OR. NEVBLK > NVECTOT-NVECPAT .OR. &
             IATJP < 1 .OR. ABS(IASPA) /= 1) CALL density_fail('Invalid mixing block')
         WRITE (*, '(5I8)') NB, NCFBLK, NEVBLK, IATJP, IASPA
         IF (JB /= NB) CALL density_fail('Mixing blocks out of order')

         IF (NEVBLK > 0) THEN

            READ (25, IOSTAT=IOS) (IVEC(NVECPAT + I),I=1,NEVBLK)
            IF (IOS /= 0) CALL density_fail('Truncated mixing state indices')
               ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
            IATJPO(NVECPAT+1:NEVBLK+NVECPAT) = IATJP
            IASPAR(NVECPAT+1:NEVBLK+NVECPAT) = IASPA

            READ (25, IOSTAT=IOS) EAV, (EVAL(NVECPAT+I),I=1,NEVBLK)

            IF (IOS /= 0) CALL density_fail('Truncated mixing energies')
!           ...Construct the true energy by adding up the average
            EVAL(NVECPAT+1:NEVBLK+NVECPAT) = EVAL(NVECPAT+1:NEVBLK+NVECPAT) + &
               EAV
!           ...For overal (all blocks) average energy
            EAVSUM = EAVSUM + EAV*NCFBLK
            NEAVSUM = NEAVSUM + NCFBLK

            READ (25, IOSTAT=IOS) ((EVEC(NVECSIZPAT+NCFPAT+I+(J-1)*NCFTOT),I=1,NCFBLK),J=1,&
               NEVBLK)
            IF (IOS /= 0) CALL density_fail('Truncated mixing eigenvectors')
         ENDIF

         NVECPAT = NVECPAT + NEVBLK
         NCFPAT = NCFPAT + NCFBLK
         NVECSIZPAT = NVECSIZPAT + NEVBLK*NCFTOT

      END DO

!     ...Here eav is the average energy of the blocks where at least
!        one eigenstate is calculated. It is not the averge of the
!        total Hamiltonian.

      IF (NCFPAT /= NCFTOT .OR. NVECPAT /= NVECTOT) &
         CALL density_fail('Mixing block totals do not match the header')
      IF (.NOT.ALL(ieee_is_finite(EVEC))) CALL density_fail('Non-finite mixing coefficients')
      DO I=1,NVECTOT
         IF (ABS(SUM(EVEC((I-1)*NCFTOT+1:I*NCFTOT)**2)-1.0D0) > 1.0D-6) &
            CALL density_fail('Mixing eigenvector is not normalized')
      END DO
      EAV = EAVSUM/NEAVSUM

      IF (NCFTOT /= NEAVSUM) WRITE (6, *) &
         'Not all blocks are diagonalized --- Average E ', 'not correct'

!     ...Substrct the overal average energy
      EVAL(:NVECTOT) = EVAL(:NVECTOT) - EAV

      CLOSE(25)

      NCF = NCFTOT
      NVEC = NVECTOT

      RETURN
      END SUBROUTINE GETMIXBLOCK
