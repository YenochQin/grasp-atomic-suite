!***********************************************************************
!                                                                      *
      SUBROUTINE SETOUT_GJ(NAME, NCI, GJC_DIAG, DGJC_DIAG)
!                                                                      *
!   Writes diagonal Landé-factor data to <name>.gj or <name>.cgj.     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE DEF_C,           ONLY: CVAC
      USE jlabl_C,         ONLY: LABJ=>JLBR, LABP=>JLBP
      USE nsmdat_C,        ONLY: SQN, DMOMNM, QMOMB,                  &
                                 HFSI=>SQN, HFSD=>DMOMNM, HFSQ=>QMOMB
      USE prnt_C,          ONLY: NVEC, IVEC
      USE syma_C,          ONLY: IATJPO, IASPAR
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=24), INTENT(IN) :: NAME
      INTEGER, INTENT(IN) :: NCI
      REAL(DOUBLE), INTENT(IN) :: GJC_DIAG(NVEC), DGJC_DIAG(NVEC)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IERR, JJ, K
      REAL(DOUBLE) :: FJ, GJA1, GJ, DGJ
      CHARACTER :: FILNAM*256, FORM*11, STATUS*3
!-----------------------------------------------
!
      K = INDEX(NAME,' ')
      IF (NCI == 0) THEN
         FILNAM = NAME(1:K-1)//'.cgj'
      ELSE
         FILNAM = NAME(1:K-1)//'.gj'
      ENDIF
!
      FORM = 'FORMATTED'
      STATUS = 'NEW'
      CALL OPENFL(29, FILNAM, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (6, *) 'Error when opening', FILNAM
         STOP
      ENDIF
!
      WRITE (29, 302) HFSI
      WRITE (29, 303) HFSD
      WRITE (29, 304) HFSQ
      WRITE (29, 402)
!
      DO I = 1, NVEC
         JJ = IATJPO(I)
         IF (JJ <= 1) CYCLE
         FJ = 0.5D00*DBLE(JJ - 1)
         GJA1 = SQRT(1.0D00/(FJ*(FJ + 1.0D00)))
         GJ = CVAC*GJA1*GJC_DIAG(I)
         DGJ = 0.001160D0*GJA1*DGJC_DIAG(I)
!
         WRITE (29, 403) IVEC(I), LABJ(JJ), LABP((IASPAR(I)+3)/2), GJ, DGJ, &
            GJ + DGJ
      END DO
!
      CLOSE(29)
      RETURN
!
  302 FORMAT('Nuclear spin                        ',1P,D22.15,' au')
  303 FORMAT('Nuclear magnetic dipole moment      ',1P,D22.15,' n.m.')
  304 FORMAT('Nuclear electric quadrupole moment  ',1P,D22.15,' barns')
  402 FORMAT(/,/,' Interaction constants:'/,/,' Level1  J Parity ',10X,'g_J', &
         14X,'delta g_J',11X,'total g_J'/)
  403 FORMAT(1X,1I3,5X,2A4,1P,3D20.10)
      END SUBROUTINE SETOUT_GJ
