!***********************************************************************
!                                                                      *
      SUBROUTINE DENSNEW_MPI(DOIT,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
!                                                                      *
!   MPI version of DENSNEW. Each rank handles a strided subset of the  *
!   outer IC loop and can optionally write a per-rank angular-data      *
!   fragment that the root process later merges into name.IOB.          *
!                                                                      *
!***********************************************************************
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW, NNNP
      USE debug_C
      USE decide_C
      USE DEF_C
      USE eigv_C
      USE foparm_C
      USE grid_C
      USE JLABL_C
      USE npar_C
      USE orb_C
      USE blk_C
      USE prnt_C
      USE TEILST_C
      USE BUFFER_C
      USE ris_C
      USE syma_C
      USE mpi_C,            ONLY: MYID, NPROCS
      USE prnt_C,           ONLY : NVEC
      USE alcbuf_I
      USE convrt_I
      USE itjpo_I
      USE onescalar_I
      IMPLICIT NONE
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2, &
                                                        DINT3, DINT4, &
                                                        DINT5, DINT6, &
                                                        DINT7
      INTEGER, INTENT(IN) :: DOIT
      INTEGER, PARAMETER :: KEY = KEYORB
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL_S
      REAL(DOUBLE) :: ELEMNT1, ELEMNT2, ELEMNT3, ELEMNT4, ELEMNT5,     &
                      ELEMNT6, ELEMNT7, CONTRI1, CONTRI2, CONTRI3,     &
                      CONTRI4, CONTRI5, CONTRI6, CONTRI7
      CHARACTER :: CNUM*11, CK*2
      INTEGER, DIMENSION(NNNW) :: IA_S
      INTEGER :: KA, IOPAR, INCOR, IC, LCNUM, ITJPOC, IR, IA, IB, I, J
      INTEGER :: LOC, NCONTR, LAB

      KA = 0
      IOPAR = 1
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
               ELEMNT1 = 0.0D00
               ELEMNT2 = 0.0D00
               ELEMNT3 = 0.0D00
               ELEMNT4 = 0.0D00
               ELEMNT5 = 0.0D00
               ELEMNT6 = 0.0D00
               ELEMNT7 = 0.0D00

               CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
               IF (IA .NE. 0) THEN
                  IF (IA .EQ. IB) THEN
                     NCONTR = 0
                     DO IA = 1,NW
                        IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                           NCONTR = NCONTR + 1
                           TSHELL_S(NCONTR) = TSHELL(IA)
                           IA_S(NCONTR) = IA
                           ELEMNT1 = ELEMNT1 + DINT1(IA,IA)*TSHELL(IA)
                           ELEMNT2 = ELEMNT2 + DINT2(IA,IA)*TSHELL(IA)
                           ELEMNT3 = ELEMNT3 + DINT3(IA,IA)*TSHELL(IA)
                           ELEMNT4 = ELEMNT4 + DINT4(IA,IA)*TSHELL(IA)
                           ELEMNT5 = ELEMNT5 + DINT5(IA,IA)*TSHELL(IA)
                           ELEMNT6 = ELEMNT6 + DINT6(IA,IA)*TSHELL(IA)
                           ELEMNT7 = ELEMNT7 + DINT7(IA,IA)*TSHELL(IA)
                        ENDIF
                     END DO
                     IF (DOIT.EQ.1) WRITE(50) IC,IR,NCONTR
                     IF (DOIT.EQ.1) THEN
                        DO I = 1,NCONTR
                           LAB = IA_S(I)*(KEY + 1)
                           WRITE(50) TSHELL_S(I),LAB
                        END DO
                     ENDIF
                  ELSE
                     IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                        IF (NAK(IA).EQ.NAK(IB)) THEN
                           IF (DOIT.EQ.1) WRITE(50) IC,IR,1
                           IF (DOIT.EQ.1) THEN
                              LAB = IA*KEY + IB
                              WRITE(50) TSHELL(1),LAB
                           ENDIF
                           ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL(1)
                           ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL(1)
                           ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL(1)
                           ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL(1)
                           ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL(1)
                           ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL(1)
                           ELEMNT7 = ELEMNT7 + DINT7(IA,IB)*TSHELL(1)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
               DO J = 1,NVEC
                  LOC = (J-1)*NCF
                  CONTRI1 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1
                  CONTRI2 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT2
                  CONTRI3 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT3
                  CONTRI4 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT4
                  CONTRI5 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT5
                  CONTRI6 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT6
                  CONTRI7 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT7
                  IF (IR.NE.IC) THEN
                     CONTRI1 = 2.0D00 * CONTRI1
                     CONTRI2 = 2.0D00 * CONTRI2
                     CONTRI3 = 2.0D00 * CONTRI3
                     CONTRI4 = 2.0D00 * CONTRI4
                     CONTRI5 = 2.0D00 * CONTRI5
                     CONTRI6 = 2.0D00 * CONTRI6
                     CONTRI7 = 2.0D00 * CONTRI7
                  ENDIF
                  DENS1(J) = DENS1(J) + CONTRI1
                  DENS2(J) = DENS2(J) + CONTRI2
                  DENS3(J) = DENS3(J) + CONTRI3
                  DENS4(J) = DENS4(J) + CONTRI4
                  DENS5(J) = DENS5(J) + CONTRI5
                  DENS6(J) = DENS6(J) + CONTRI6
                  DENS7(J) = DENS7(J) + CONTRI7
               END DO
            ENDIF
         END DO
      END DO

      CALL ALCBUF (3)
      RETURN
      END SUBROUTINE DENSNEW_MPI
