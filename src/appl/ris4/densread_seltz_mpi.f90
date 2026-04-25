!***********************************************************************
!                                                                      *
      SUBROUTINE DENSREAD_SELTZ_MPI(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,&
         DINT7,DINT1VEC,DENS1VEC,NRNUC)
!                                                                      *
!   MPI version of DENSREAD_SELTZ. Each rank scans the full angular    *
!   file but only accumulates records whose IC belongs to its strided   *
!   subset.                                                             *
!                                                                      *
!***********************************************************************
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW, NNNP
      USE debug_C,          ONLY: CUTOFF
      USE prnt_C
      USE ris_C
      USE orb_C
      USE eigv_C
      USE mpi_C,            ONLY: MYID, NPROCS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NRNUC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2,  &
                                          DINT3, DINT4, DINT5, DINT6,  &
                                          DINT7
      REAL(DOUBLE), DIMENSION(NVEC,NRNUC), INTENT(OUT) :: DENS1VEC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW,NRNUC), INTENT(IN) :: DINT1VEC
      INTEGER, PARAMETER :: KEY = KEYORB
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL_R
      REAL(DOUBLE), DIMENSION(NRNUC) :: CONTRI1VEC, ELEMNT1VEC
      REAL(DOUBLE) :: ELEMNT1, ELEMNT2, ELEMNT3, ELEMNT4, ELEMNT5
      REAL(DOUBLE) :: ELEMNT6, ELEMNT7
      REAL(DOUBLE) :: CONTRI1, CONTRI2, CONTRI3, CONTRI4, CONTRI5
      REAL(DOUBLE) :: CONTRI6, CONTRI7
      REAL(DOUBLE) :: EVPAIR, PAIR_SCALE
      INTEGER :: IOS, IA, IB, IC, IR, I, J, L, LOC, LAB, NCOUNT
      LOGICAL :: OWNREC

      DENS1VEC(:,:) = 0.0D00
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
         ELEMNT1VEC(:) = 0.0D00

         DO I= 1,NCOUNT
            READ(50) TSHELL_R(I),LAB
            IF (.NOT. OWNREC) CYCLE
            IA = LAB/KEY
            IB = MOD(LAB,KEY)
            ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL_R(I)
            DO L = 2,NRNUC
               ELEMNT1VEC(L) = ELEMNT1VEC(L) + DINT1VEC(IA,IB,L)*TSHELL_R(I)
            END DO
            ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL_R(I)
            ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL_R(I)
            ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL_R(I)
            ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL_R(I)
            ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL_R(I)
            ELEMNT7 = ELEMNT7 + DINT7(IA,IB)*TSHELL_R(I)
         ENDDO

         IF (OWNREC) THEN
            IF (ABS(ELEMNT1) .LE. CUTOFF .AND. &
                ABS(ELEMNT2) .LE. CUTOFF .AND. &
                ABS(ELEMNT3) .LE. CUTOFF .AND. &
                ABS(ELEMNT4) .LE. CUTOFF .AND. &
                ABS(ELEMNT5) .LE. CUTOFF .AND. &
                ABS(ELEMNT6) .LE. CUTOFF .AND. &
                ABS(ELEMNT7) .LE. CUTOFF .AND. &
                ALL(ABS(ELEMNT1VEC(:)) .LE. CUTOFF)) GOTO 16
            IF (IR .EQ. IC) THEN
               PAIR_SCALE = 1.0D00
            ELSE
               PAIR_SCALE = 2.0D00
            ENDIF
            DO J = 1,NVEC
               LOC = (J-1)*NCF
               EVPAIR = PAIR_SCALE * EVEC(IC+LOC)*EVEC(IR+LOC)
               CONTRI1 = EVPAIR*ELEMNT1
               CONTRI1VEC(:) = EVPAIR*ELEMNT1VEC(:)
               CONTRI2 = EVPAIR*ELEMNT2
               CONTRI3 = EVPAIR*ELEMNT3
               CONTRI4 = EVPAIR*ELEMNT4
               CONTRI5 = EVPAIR*ELEMNT5
               CONTRI6 = EVPAIR*ELEMNT6
               CONTRI7 = EVPAIR*ELEMNT7
               DENS1(J) = DENS1(J) + CONTRI1
               DO L = 2,NRNUC
                  DENS1VEC(J,L) = DENS1VEC(J,L) + CONTRI1VEC(L)
               END DO
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
      END SUBROUTINE DENSREAD_SELTZ_MPI
