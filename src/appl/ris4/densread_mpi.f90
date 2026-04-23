!***********************************************************************
!                                                                      *
      SUBROUTINE DENSREAD_MPI(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
!                                                                      *
!   MPI version of DENSREAD. Each rank scans the full angular file but *
!   only accumulates records whose IC belongs to its strided subset.    *
!                                                                      *
!***********************************************************************
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW
      USE prnt_C
      USE ris_C
      USE orb_C
      USE eigv_C
      USE mpi_C,            ONLY: MYID, NPROCS
      IMPLICIT NONE
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2,  &
                                          DINT3, DINT4, DINT5, DINT6,  &
                                          DINT7
      INTEGER, PARAMETER :: KEY = KEYORB
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL_R
      REAL(DOUBLE) :: ELEMNT1, ELEMNT2, ELEMNT3, ELEMNT4, ELEMNT5
      REAL(DOUBLE) :: ELEMNT6, ELEMNT7
      REAL(DOUBLE) :: CONTRI1, CONTRI2, CONTRI3, CONTRI4, CONTRI5
      REAL(DOUBLE) :: CONTRI6, CONTRI7
      INTEGER :: IOS, IA, IB, IC, IR, I, J, LOC, LAB, NCOUNT
      LOGICAL :: OWNREC

      REWIND (50)
   16 READ (50,IOSTAT = IOS) IC,IR,NCOUNT

      IF (IOS .EQ. 0) THEN
         OWNREC = MOD(IC-1, NPROCS) .EQ. MYID
         ELEMNT1 = 0.0D00
         ELEMNT2 = 0.0D00
         ELEMNT3 = 0.0D00
         ELEMNT4 = 0.0D00
         ELEMNT5 = 0.0D00
         ELEMNT6 = 0.0D00
         ELEMNT7 = 0.0D00

         DO I= 1,NCOUNT
            READ(50) TSHELL_R(I),LAB
            IF (.NOT. OWNREC) CYCLE
            IA = LAB/KEY
            IB = MOD(LAB,KEY)
            ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL_R(I)
            ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL_R(I)
            ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL_R(I)
            ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL_R(I)
            ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL_R(I)
            ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL_R(I)
            ELEMNT7 = ELEMNT7 + DINT7(IA,IB)*TSHELL_R(I)
         ENDDO

         IF (OWNREC) THEN
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

         GOTO 16
      ENDIF

      RETURN
      END SUBROUTINE DENSREAD_MPI
