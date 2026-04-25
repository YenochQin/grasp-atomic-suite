!***********************************************************************
!                                                                      *
      SUBROUTINE SMSREAD_MPI(VINT,VINT2)
!                                                                      *
!   MPI version of SMSREAD. Each rank scans the full angular file but  *
!   only accumulates records whose IC belongs to its strided subset.    *
!                                                                      *
!***********************************************************************
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW
      USE prnt_C
      USE ris_C
      USE orb_C
      USE eigv_C
      USE BUFFER_C
      USE debug_C,          ONLY: CUTOFF
      USE mpi_C,            ONLY: MYID, NPROCS
      USE alcbuf_I
      IMPLICIT NONE
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: VINT, VINT2
      INTEGER, PARAMETER :: KEY = KEYORB
      REAL(DOUBLE), DIMENSION(NVEC) :: EVPAIR
      REAL(DOUBLE) :: CONTRI, CONTRIK1, COEFFSMS, RADIAL_K1, RADIAL_SMS,&
         PAIR_SCALE, TERM_K1, TERM_SMS
      INTEGER :: IC, IR, J, IIA, IIB, IIC, IID, LOC, LAB, IOS
      LOGICAL :: OWNREC

      CALL ALCBUF (1)
      REWIND(51)
   16 READ (51,IOSTAT = IOS) IC,IR
      IF (IOS .EQ. 0) THEN
         READ(51) COEFFSMS,LAB
         OWNREC = MOD(IC-1, NPROCS) .EQ. MYID
         IF (OWNREC) THEN
            IID = MOD (LAB, KEY)
            LAB = LAB/KEY
            IIB  = MOD (LAB, KEY)
            LAB = LAB/KEY
            IIC = MOD (LAB, KEY)
            IIA = LAB/KEY
            RADIAL_K1 = VINT(IIA,IIC) * VINT(IIB,IID)
            RADIAL_SMS = (VINT2(IIA,IIC) * VINT(IIB,IID)                 &
               + VINT2(IIB,IID) * VINT(IIA,IIC))/2.0D00
            IF (RADIAL_K1 == 0.0D00 .AND. RADIAL_SMS == 0.0D00) GOTO 16
            TERM_K1 = -COEFFSMS * RADIAL_K1
            TERM_SMS = -COEFFSMS * RADIAL_SMS
            IF (IR .EQ. IC) THEN
               PAIR_SCALE = 1.0D00
            ELSE
               PAIR_SCALE = 2.0D00
            ENDIF
            DO J = 1,NVEC
               LOC = (J-1)*NCF
               EVPAIR(J) = PAIR_SCALE * EVEC(IC+LOC) * EVEC(IR+LOC)
            END DO
            DO J = 1,NVEC
               CONTRIK1 = EVPAIR(J) * TERM_K1
               CONTRI = EVPAIR(J) * TERM_SMS
               SMSC1(J) = SMSC1(J) + CONTRIK1
               SMSC2(J) = SMSC2(J) + CONTRI
            ENDDO
         ENDIF
         GOTO 16
      ENDIF
      CALL ALCBUF (3)
      RETURN
      END SUBROUTINE SMSREAD_MPI
