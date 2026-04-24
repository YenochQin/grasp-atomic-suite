!*******************************************************************
!                                                                  *
      MODULE oneparticlejj2_stats_C
!                                                                  *
!*******************************************************************
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def, ONLY: NNNW
      USE m_C, ONLY: NPEEL, JLIST, JJC1, JJC2, JJQ1, JJQ2
      INTEGER, PARAMETER :: OPJJ2_CACHE_SIZE = 256
      INTEGER :: OPJJ2_CALLS = 0
      INTEGER :: OPJJ2_FAIL_RECOP00 = 0
      INTEGER :: OPJJ2_FAIL_RECOP2_PRE = 0
      INTEGER :: OPJJ2_FAIL_IK1 = 0
      INTEGER :: OPJJ2_FAIL_IK2 = 0
      INTEGER :: OPJJ2_FAIL_C0T5S_1 = 0
      INTEGER :: OPJJ2_FAIL_C0T5S_2 = 0
      INTEGER :: OPJJ2_FAIL_RMEAJJ_1 = 0
      INTEGER :: OPJJ2_FAIL_RMEAJJ_2 = 0
      INTEGER :: OPJJ2_SUCCESS = 0
      INTEGER :: OPJJ2_CACHE_HIT = 0
      INTEGER :: OPJJ2_CACHE_MISS = 0
      LOGICAL :: OPJJ2_CACHE_VALID(OPJJ2_CACHE_SIZE) = .FALSE.
      INTEGER :: OPJJ2_CACHE_STAGE(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_NS(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_JA(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_JB(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_KA(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_K1(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_K2(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_NPEEL(OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_JLIST(NNNW,OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_JJC1(NNNW,OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_JJC2(NNNW,OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_JJQ1(3,NNNW,OPJJ2_CACHE_SIZE) = 0
      INTEGER :: OPJJ2_CACHE_JJQ2(3,NNNW,OPJJ2_CACHE_SIZE) = 0
      REAL(DOUBLE) :: OPJJ2_CACHE_RECC(OPJJ2_CACHE_SIZE) = 0.0D00
      CONTAINS
      SUBROUTINE RESET_ONEPARTICLEJJ2_STATS()
      OPJJ2_CALLS = 0
      OPJJ2_FAIL_RECOP00 = 0
      OPJJ2_FAIL_RECOP2_PRE = 0
      OPJJ2_FAIL_IK1 = 0
      OPJJ2_FAIL_IK2 = 0
      OPJJ2_FAIL_C0T5S_1 = 0
      OPJJ2_FAIL_C0T5S_2 = 0
      OPJJ2_FAIL_RMEAJJ_1 = 0
      OPJJ2_FAIL_RMEAJJ_2 = 0
      OPJJ2_SUCCESS = 0
      OPJJ2_CACHE_HIT = 0
      OPJJ2_CACHE_MISS = 0
      OPJJ2_CACHE_VALID(:) = .FALSE.
      END SUBROUTINE RESET_ONEPARTICLEJJ2_STATS
      INTEGER FUNCTION ONEPARTICLEJJ2_CACHE_SLOT(NS,JA,JB,KA,K1,K2)
      INTEGER, INTENT(IN) :: NS, JA, JB, KA, K1, K2
      INTEGER :: HASHV, I
      HASHV = MOD(NS + 31*JA + 37*JB + 41*KA + 43*K1 + 47*K2 + 53*NPEEL, &
         OPJJ2_CACHE_SIZE)
      IF (HASHV < 0) HASHV = HASHV + OPJJ2_CACHE_SIZE
      DO I = 1, NPEEL
         HASHV = MOD(HASHV*131 + JLIST(I)*7 + JJQ1(1,JLIST(I))*11 + &
            JJQ1(2,JLIST(I))*13 + JJQ1(3,JLIST(I))*17 + JJQ2(1,JLIST(I))*19 + &
            JJQ2(2,JLIST(I))*23 + JJQ2(3,JLIST(I))*29, OPJJ2_CACHE_SIZE)
         IF (HASHV < 0) HASHV = HASHV + OPJJ2_CACHE_SIZE
      END DO
      IF (NPEEL > 1) THEN
         DO I = 1, NPEEL - 1
            HASHV = MOD(HASHV*131 + JJC1(I)*31 + JJC2(I)*37, OPJJ2_CACHE_SIZE)
            IF (HASHV < 0) HASHV = HASHV + OPJJ2_CACHE_SIZE
         END DO
      ENDIF
      ONEPARTICLEJJ2_CACHE_SLOT = HASHV + 1
      END FUNCTION ONEPARTICLEJJ2_CACHE_SLOT
      LOGICAL FUNCTION ONEPARTICLEJJ2_CACHE_MATCH(SLOT,NS,JA,JB,KA,K1,K2)
      INTEGER, INTENT(IN) :: SLOT, NS, JA, JB, KA, K1, K2
      INTEGER :: I, IJ
      ONEPARTICLEJJ2_CACHE_MATCH = .FALSE.
      IF (.NOT. OPJJ2_CACHE_VALID(SLOT)) RETURN
      IF (OPJJ2_CACHE_NS(SLOT) /= NS) RETURN
      IF (OPJJ2_CACHE_JA(SLOT) /= JA) RETURN
      IF (OPJJ2_CACHE_JB(SLOT) /= JB) RETURN
      IF (OPJJ2_CACHE_KA(SLOT) /= KA) RETURN
      IF (OPJJ2_CACHE_K1(SLOT) /= K1) RETURN
      IF (OPJJ2_CACHE_K2(SLOT) /= K2) RETURN
      IF (OPJJ2_CACHE_NPEEL(SLOT) /= NPEEL) RETURN
      DO I = 1, NPEEL
         IJ = JLIST(I)
         IF (OPJJ2_CACHE_JLIST(I,SLOT) /= IJ) RETURN
         IF (ANY(OPJJ2_CACHE_JJQ1(:,I,SLOT) /= JJQ1(:,IJ))) RETURN
         IF (ANY(OPJJ2_CACHE_JJQ2(:,I,SLOT) /= JJQ2(:,IJ))) RETURN
      END DO
      IF (NPEEL > 1) THEN
         DO I = 1, NPEEL - 1
            IF (OPJJ2_CACHE_JJC1(I,SLOT) /= JJC1(I)) RETURN
            IF (OPJJ2_CACHE_JJC2(I,SLOT) /= JJC2(I)) RETURN
         END DO
      ENDIF
      ONEPARTICLEJJ2_CACHE_MATCH = .TRUE.
      END FUNCTION ONEPARTICLEJJ2_CACHE_MATCH
      SUBROUTINE STORE_ONEPARTICLEJJ2_CACHE(SLOT,NS,JA,JB,KA,K1,K2,STAGE,RECC)
      INTEGER, INTENT(IN) :: SLOT, NS, JA, JB, KA, K1, K2, STAGE
      REAL(DOUBLE), INTENT(IN) :: RECC
      INTEGER :: I, IJ
      OPJJ2_CACHE_VALID(SLOT) = .TRUE.
      OPJJ2_CACHE_STAGE(SLOT) = STAGE
      OPJJ2_CACHE_NS(SLOT) = NS
      OPJJ2_CACHE_JA(SLOT) = JA
      OPJJ2_CACHE_JB(SLOT) = JB
      OPJJ2_CACHE_KA(SLOT) = KA
      OPJJ2_CACHE_K1(SLOT) = K1
      OPJJ2_CACHE_K2(SLOT) = K2
      OPJJ2_CACHE_NPEEL(SLOT) = NPEEL
      OPJJ2_CACHE_JLIST(:,SLOT) = 0
      OPJJ2_CACHE_JJC1(:,SLOT) = 0
      OPJJ2_CACHE_JJC2(:,SLOT) = 0
      OPJJ2_CACHE_JJQ1(:,:,SLOT) = 0
      OPJJ2_CACHE_JJQ2(:,:,SLOT) = 0
      DO I = 1, NPEEL
         IJ = JLIST(I)
         OPJJ2_CACHE_JLIST(I,SLOT) = IJ
         OPJJ2_CACHE_JJQ1(:,I,SLOT) = JJQ1(:,IJ)
         OPJJ2_CACHE_JJQ2(:,I,SLOT) = JJQ2(:,IJ)
      END DO
      IF (NPEEL > 1) THEN
         OPJJ2_CACHE_JJC1(1:NPEEL-1,SLOT) = JJC1(1:NPEEL-1)
         OPJJ2_CACHE_JJC2(1:NPEEL-1,SLOT) = JJC2(1:NPEEL-1)
      ENDIF
      OPJJ2_CACHE_RECC(SLOT) = RECC
      END SUBROUTINE STORE_ONEPARTICLEJJ2_CACHE
      END MODULE oneparticlejj2_stats_C
!*******************************************************************
!                                                                  *
      SUBROUTINE ONEPARTICLEJJ2(NS,KA,JA,JB,COEFF)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1 +- 1   *
!                                                  N'2 = N2 -+ 1   *
!                                                                  *
!      SUBROUTINE CALLED:                                          *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, TENTH, HALF, EPS
      USE m_C,             ONLY: NQ1, JLIST
      USE orb_C,           ONLY: NAK
      USE oneparticlejj2_stats_C
      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recop00_I
      USE recop2_I
      USE c0t5s_I
      USE rmeajj_I
      USE wj1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: NS,KA,JA,JB
      REAL(DOUBLE), INTENT(OUT) :: COEFF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: I,IA,IB,IAT,IJ,IFAZ,IQMM1,IQMM2,JIBKS1,JIB, &
                      KS1,KS2,CACHE_SLOT,CACHE_STAGE
      REAL(DOUBLE) :: REC,A1,A2,A3,S1,S2,QM1,QM2
!-----------------------------------------------
!
!     THE CASE 12   + -
!
      OPJJ2_CALLS = OPJJ2_CALLS + 1
      COEFF=ZERO
      IA=MIN0(JA,JB)
      IB=MAX0(JA,JB)
      IJ=JLIST(IA)
      KS1=(IABS(NAK(IJ))*2)-1
      IJ=JLIST(IB)
      KS2=(IABS(NAK(IJ))*2)-1
      CACHE_SLOT = ONEPARTICLEJJ2_CACHE_SLOT(NS,JA,JB,KA,KS1,KS2)
      IF (ONEPARTICLEJJ2_CACHE_MATCH(CACHE_SLOT,NS,JA,JB,KA,KS1,KS2)) THEN
         OPJJ2_CACHE_HIT = OPJJ2_CACHE_HIT + 1
         CACHE_STAGE = OPJJ2_CACHE_STAGE(CACHE_SLOT)
         IF (CACHE_STAGE == 0) THEN
            OPJJ2_FAIL_RECOP00 = OPJJ2_FAIL_RECOP00 + 1
            RETURN
         ELSE IF (CACHE_STAGE == 1) THEN
            OPJJ2_FAIL_RECOP2_PRE = OPJJ2_FAIL_RECOP2_PRE + 1
            RETURN
         ELSE
            REC = OPJJ2_CACHE_RECC(CACHE_SLOT)
            GOTO 10
         ENDIF
      ENDIF
      OPJJ2_CACHE_MISS = OPJJ2_CACHE_MISS + 1
      CALL RECOP00(NS,IA,IB,KA,IAT)
      IF(IAT == 0) THEN
         OPJJ2_FAIL_RECOP00 = OPJJ2_FAIL_RECOP00 + 1
         CALL STORE_ONEPARTICLEJJ2_CACHE(CACHE_SLOT,NS,JA,JB,KA,KS1,KS2,0,ZERO)
         RETURN
      ENDIF
      CALL RECOP2(NS,IA,IB,KS1,KS2,KA,0,IAT,REC)
      IF(IAT == 0) THEN
         OPJJ2_FAIL_RECOP2_PRE = OPJJ2_FAIL_RECOP2_PRE + 1
         CALL STORE_ONEPARTICLEJJ2_CACHE(CACHE_SLOT,NS,JA,JB,KA,KS1,KS2,1,ZERO)
         RETURN
      ENDIF
      CALL STORE_ONEPARTICLEJJ2_CACHE(CACHE_SLOT,NS,JA,JB,KA,KS1,KS2,2,REC)
   10 CONTINUE
      CALL PERKO2(JA,JB,JA,JA,2)
      QM1=HALF
      QM2=-HALF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4) /= (ID1(4)+IQMM1)) THEN
         OPJJ2_FAIL_IK1 = OPJJ2_FAIL_IK1 + 1
         RETURN
      ENDIF
      IQMM2=QM2+QM2+TENTH*QM2
      IF(IK2(4) /= (ID2(4)+IQMM2)) THEN
         OPJJ2_FAIL_IK2 = OPJJ2_FAIL_IK2 + 1
         RETURN
      ENDIF
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A2)
      IF(DABS(A2) < EPS) THEN
         OPJJ2_FAIL_C0T5S_1 = OPJJ2_FAIL_C0T5S_1 + 1
         RETURN
      ENDIF
      CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),A3)
      IF(DABS(A3) < EPS) THEN
         OPJJ2_FAIL_C0T5S_2 = OPJJ2_FAIL_C0T5S_2 + 1
         RETURN
      ENDIF
      CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S1)
      IF(DABS(S1) < EPS) THEN
         OPJJ2_FAIL_RMEAJJ_1 = OPJJ2_FAIL_RMEAJJ_1 + 1
         RETURN
      ENDIF
      CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S2)
      IF(DABS(S2) < EPS) THEN
         OPJJ2_FAIL_RMEAJJ_2 = OPJJ2_FAIL_RMEAJJ_2 + 1
         RETURN
      ENDIF
      A1=S1*S2*A2*A3
      CALL RECOP2(NS,IA,IB,KS1,KS2,KA,1,IAT,REC)
      COEFF=A1*REC/DSQRT(DBLE((2*KA+1)*(IK1(7)+1)*(IK2(7)+1)))
      JIB=IB-1
      IFAZ=0
      DO  I=IA,JIB
        IJ=JLIST(I)
        IFAZ=IFAZ+NQ1(IJ)
      END DO
      IFAZ=IFAZ+1
      IF(MOD(IFAZ,2) /= 0)COEFF=-COEFF
      IF(IA.NE.JA) THEN
        IFAZ = IK1(3)+IK2(3)-2*KA+2
        IF(MOD(IFAZ,4) /= 0)COEFF=-COEFF
      ENDIF
      COEFF=-COEFF*SQRT(DBLE(ID1(3)+1))
      OPJJ2_SUCCESS = OPJJ2_SUCCESS + 1
      RETURN
      END SUBROUTINE ONEPARTICLEJJ2
