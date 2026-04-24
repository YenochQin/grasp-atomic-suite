!***********************************************************************
!                                                                      *
      SUBROUTINE SMSNEW_MPI(DOIT,VINT,VINT2)
!                                                                      *
!   MPI version of SMSNEW. Each rank handles a strided subset of the   *
!   outer IC loop and can optionally write a per-rank angular-data      *
!   fragment that the root process later merges into name.ITB.          *
!                                                                      *
!***********************************************************************
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW
      USE debug_C
      USE prnt_C
      USE orb_C
      USE BUFFER_C
      USE eigv_C
      USE ris_C
      USE mpi_C,            ONLY: MYID, NPROCS
      USE alcbuf_I
      USE itjpo_I
      IMPLICIT NONE
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: VINT, VINT2
      INTEGER, INTENT(IN) :: DOIT
      EXTERNAL :: CORD
      INTEGER, PARAMETER :: KEY = KEYORB
      CHARACTER*11 :: CNUM
      REAL(DOUBLE) :: CONTRI, CONTRIK1, VCOEFF
      INTEGER :: J,K,LOC,IC,IR,ITJPOC,INCOR,LCNUM

      INCOR = 1
      CALL ALCBUF (1)

      DO IC = MYID + 1, NCF, NPROCS
         IF (MYID == 0 .AND. MOD(IC,100) .EQ. 1) THEN
            CALL CONVRT (IC,CNUM,LCNUM)
            PRINT *, 'Column '//CNUM(1:LCNUM)//' complete on root;'
         ENDIF

         ITJPOC = ITJPO (IC)
         DO IR = IC,NCF
            IF (ITJPO(IR) .EQ. ITJPOC) THEN
               NVCOEF = 0
               CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
               DO K = 1,NVCOEF
                  VCOEFF = COEFF(K)
                  IF (ABS (VCOEFF) .GT. CUTOFF) THEN
                     IF (LABEL(5,K) .EQ. 1) THEN
                        IF(DOIT.EQ.1) WRITE(51) IC,IR
                        IF(DOIT.EQ.1) THEN
                           WRITE(51) VCOEFF, ((LABEL(1,K)*KEY + LABEL(3,K))*KEY+LABEL(2,K))*KEY+LABEL(4,K)
                        ENDIF
                        DO J = 1,NVEC
                           LOC = (J-1)*NCF
                           CONTRIK1 = - EVEC(IC+LOC)*EVEC(IR+LOC)        &
                              * VCOEFF                                   &
                              * VINT (LABEL(1,K),LABEL(3,K))             &
                              * VINT (LABEL(2,K),LABEL(4,K))
                           CONTRI = - EVEC(IC+LOC)*EVEC(IR+LOC)          &
                              * VCOEFF                                   &
                              * ( VINT2(LABEL(1,K),LABEL(3,K))           &
                              * VINT (LABEL(2,K),LABEL(4,K))             &
                              + VINT2(LABEL(2,K),LABEL(4,K))             &
                              * VINT (LABEL(1,K),LABEL(3,K)) )/2.0D00
                           IF (IR.NE.IC) THEN
                              CONTRI = 2.0D00 * CONTRI
                              CONTRIK1 = 2.0D00 * CONTRIK1
                           ENDIF
                           SMSC1(J) = SMSC1(J) + CONTRIK1
                           SMSC2(J) = SMSC2(J) + CONTRI
                        END DO
                     ENDIF
                  ENDIF
               END DO
            ENDIF
         END DO
      END DO

      CALL ALCBUF (3)
      RETURN
      END SUBROUTINE SMSNEW_MPI
