!***********************************************************************
!                                                                      *
      SUBROUTINE RIS_CAL_MPI (NAME)
!                                                                      *
!   MPI version of RIS_CAL. Each rank accumulates a strided subset of  *
!   configuration-pair contributions and the partial sums are reduced   *
!   before the root rank writes the final output.                       *
!                                                                      *
!***********************************************************************
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE memory_man
      USE debug_C
      USE decide_c
      USE def_C
      USE grid_C
      USE npar_C
      USE prnt_C
      USE syma_C
      USE orb_C
      USE teilst_C
      USE buffer_C
      USE ris_C
      USE jlabl_C, LABJ=> JLBR, LABP=>JLBP
      USE eigv_C
      USE mpi_C
      USE rintdens_I
      USE rinti_I
      USE rint_I
      USE rintdensvec_I
      USE rinti_nms_I
      USE cre_I
      USE vinti_I
      USE rint_sms2_I
      USE rint_sms3_I
      USE angdata_I
      USE getyn_I
      USE densread_I
      USE densread_seltz_I
      USE densnew_seltz_I
      USE densread_mpi_I
      USE densread_seltz_mpi_I
      USE densnew_I
      USE densnew_mpi_I
      USE densnew_seltz_mpi_I
      USE smsread_I
      USE smsread_mpi_I
      USE smsnew_I
      USE smsnew_mpi_I
      USE edensityfit_I
      IMPLICIT NONE
      CHARACTER*24, INTENT(IN) :: NAME
      REAL(DOUBLE), DIMENSION(:,:,:), pointer :: DINT1VEC
      REAL(DOUBLE), DIMENSION(:,:), pointer   :: DENS1VEC
      REAL(DOUBLE), DIMENSION(:), pointer     :: DENSFIT
      REAL(DOUBLE), DIMENSION(:,:), pointer   :: FMAT
      REAL(DOUBLE), DIMENSION(:), pointer     :: RHO
      REAL(DOUBLE), DIMENSION(:), pointer     :: RES
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: VINT, VINT2
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: DINT1, DINT2, DINT3, &
                                            DINT4, DINT5, DINT6, &
                                            DINT7
      REAL(DOUBLE) :: AU2FM, RCRE, HzSMSu, HzNMSu, FO90
      LOGICAL :: AVAIL_TB, AVAIL_OB, YES2, SAVEYES, SERIAL_SAVE
      INTEGER :: I, J, ncount1, DOIT_OB, DOIT_TB, NRNUC, YES2_MPI
      INTEGER :: AVAIL_OB_MPI, AVAIL_TB_MPI, SERIAL_SAVE_MPI
      IF (MYID .EQ. 0) THEN
         WRITE(*,*) '-------------------------------'
         WRITE(*,*) 'RIS_CAL_MPI: Execution Begins ...'
         WRITE(*,*) '-------------------------------'
      ENDIF

      DOIT_OB = 0
      DOIT_TB = 0
      YES2_MPI = 0
      AVAIL_OB_MPI = 0
      AVAIL_TB_MPI = 0
      SERIAL_SAVE_MPI = 0

      FO90 = PARM(1)+2.d0*LOG(3.d0)*PARM(2)
      DO I=1,NNNP
         IF(R(I+1).GT.FO90.AND.R(I).LE.FO90) THEN
            NRNUC = I+1
            EXIT
         END IF
      END DO

      ALLOCATE( DENS1VEC(NVEC,NRNUC),DINT1VEC(NNNW,NNNW,NRNUC) )
      ALLOCATE( DENSFIT(NRNUC) )
      DINT1VEC(:,:,:) = 0.D0
      DENS1VEC(:,:) = 0.D0
      ALLOCATE( FMAT(6,NVEC),RHO(NVEC),RES(NVEC) )

      CALL ALLOC (SMSC1, NVEC,'SMSC1', 'RIS_CAL_MPI')
      CALL ALLOC (SMSC2, NVEC,'SMSC2', 'RIS_CAL_MPI')
      CALL ALLOC (DENS1, NVEC,'DENS1','RIS_CAL_MPI')
      CALL ALLOC (DENS2, NVEC,'DENS2','RIS_CAL_MPI')
      CALL ALLOC (DENS3, NVEC,'DENS3','RIS_CAL_MPI')
      CALL ALLOC (DENS4, NVEC,'DENS4','RIS_CAL_MPI')
      CALL ALLOC (DENS5, NVEC,'DENS5','RIS_CAL_MPI')
      CALL ALLOC (DENS6, NVEC,'DENS6','RIS_CAL_MPI')
      CALL ALLOC (DENS7, NVEC,'DENS7','RIS_CAL_MPI')

      CALL STARTTIME (ncount1, 'RIS_CAL_MPI')

      DO I = 1,NVEC
         SMSC1(I)  = 0.0D00
         SMSC2(I)  = 0.0D00
         DENS1(I) = 0.0D00
         DENS2(I) = 0.0D00
         DENS3(I) = 0.0D00
         DENS4(I) = 0.0D00
         DENS5(I) = 0.0D00
         DENS6(I) = 0.0D00
         DENS7(I) = 0.0D00
      ENDDO

      IF (MYID .EQ. 0) THEN
         WRITE(*,*) ' Compute higher order field shift electronic factors?'
         YES2 = GETYN ()
         IF (YES2) THEN
            YES2_MPI = 1
         ELSE
            YES2_MPI = 0
         ENDIF
      ENDIF
      CALL MPI_BCAST(YES2_MPI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      YES2 = YES2_MPI .EQ. 1

      DO I = 1,NW
         DO J = 1,NW
            IF (NAK(I).EQ.NAK(J)) THEN
               DINT1(I,J) = RINTDENS(I,J)
               IF(YES2) THEN
                  CALL RINTDENSVEC(I,J,DINT1VEC,NRNUC)
               END IF
               CALL RINTI_NMS(I,J,DINT2(I,J),DINT7(I,J))
               DINT3(I,J) = RINT(I,J,1)
               DINT4(I,J) = RINT(I,J,2)
               DINT5(I,J) = RINT(I,J,-1)
               DINT6(I,J) = RINT(I,J,-2)
            ELSE
               DINT1(I,J) = 0.0D00
               DINT1VEC(I,J,:) = 0.0D00
               DINT2(I,J) = 0.0D00
               DINT3(I,J) = 0.0D00
               DINT4(I,J) = 0.0D00
               DINT5(I,J) = 0.0D00
               DINT6(I,J) = 0.0D00
               DINT7(I,J) = 0.0D00
            ENDIF
         END DO
      END DO

      DO I = 1,NW
         DO J = 1,NW
            IF (I.NE.J) THEN
               RCRE = CRE(NAK(I),1,NAK(J))
               IF (DABS(RCRE) .GT. CUTOFF) THEN
                  VINT  (I,J) = VINTI(I,J)
                  VINT2(I,J) = VINT(I,J)                                  &
                     + RINT_SMS2(I,J)/RCRE                                &
                     + RINT_SMS3(I,J)
               ELSE
                  VINT (I,J) = 0.0D00
                  VINT2(I,J) = 0.0D00
               ENDIF
            ELSE
               VINT (I,J) = 0.0D00
               VINT2(I,J) = 0.0D00
            ENDIF
         END DO
      END DO

      IF (MYID .EQ. 0) THEN
         CALL ANGDATA(NAME,AVAIL_OB,1)
         IF (AVAIL_OB) CLOSE(50)
         CALL ANGDATA(NAME,AVAIL_TB,2)
         IF (AVAIL_TB) CLOSE(51)
         AVAIL_OB_MPI = MERGE(1, 0, AVAIL_OB)
         AVAIL_TB_MPI = MERGE(1, 0, AVAIL_TB)
         IF ((.NOT. AVAIL_OB) .AND. (.NOT. AVAIL_TB)) THEN
            PRINT *,' Save ang. coefficients of one- and two-body op.?'
            SAVEYES = GETYN ()
            PRINT *
            IF(SAVEYES) THEN
               DOIT_OB = 1
               DOIT_TB = 1
            ENDIF
         ELSEIF (.NOT. AVAIL_OB) THEN
            PRINT *,' Save ang. coefficients of one-body op. ?'
            SAVEYES = GETYN ()
            PRINT *
            IF(SAVEYES) DOIT_OB = 1
         ELSEIF (.NOT. AVAIL_TB) THEN
            PRINT *,' Save ang. coefficients of two-body op. ?'
            SAVEYES = GETYN ()
            PRINT *
            IF(SAVEYES) DOIT_TB = 1
         ENDIF
      ENDIF
      CALL MPI_BCAST(AVAIL_OB_MPI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      CALL MPI_BCAST(AVAIL_TB_MPI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      CALL MPI_BCAST(DOIT_OB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      CALL MPI_BCAST(DOIT_TB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      AVAIL_OB = AVAIL_OB_MPI .EQ. 1
      AVAIL_TB = AVAIL_TB_MPI .EQ. 1
      IF (MYID .EQ. 0) THEN
         ! Saving angular coefficients still has to stay on the root
         ! process because the serial file format is shared.
         SERIAL_SAVE = ((.NOT. AVAIL_OB .AND. DOIT_OB .EQ. 1) .OR.      &
                        (.NOT. AVAIL_TB .AND. DOIT_TB .EQ. 1))
         SERIAL_SAVE_MPI = MERGE(1, 0, SERIAL_SAVE)
      ENDIF
      CALL MPI_BCAST(SERIAL_SAVE_MPI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
      SERIAL_SAVE = SERIAL_SAVE_MPI .EQ. 1

      IF (SERIAL_SAVE) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
         IF (MYID .EQ. 0) THEN
            J = INDEX(NAME,' ')
            IF (AVAIL_OB) THEN
               OPEN(UNIT=50,FILE = NAME(1:J-1)//'.IOB',STATUS='OLD',FORM='UNFORMATTED')
               IF (YES2) THEN
                  CALL DENSREAD_SELTZ(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6, &
                       DINT7,DINT1VEC,DENS1VEC,NRNUC)
               ELSE
                  CALL DENSREAD(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
               ENDIF
               CLOSE(50)
            ELSE
               IF (DOIT_OB .EQ. 1) THEN
                  OPEN(UNIT=50,FILE = NAME(1:J-1)//'.IOB',STATUS='UNKNOWN',FORM='UNFORMATTED')
               ENDIF
               IF (YES2) THEN
                  CALL DENSNEW_SELTZ(DOIT_OB,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6, &
                       DINT7,DINT1VEC,DENS1VEC,NRNUC)
               ELSE
                  CALL DENSNEW(DOIT_OB,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
               ENDIF
            ENDIF

            IF (AVAIL_TB) THEN
               OPEN(UNIT=51,FILE = NAME(1:J-1)//'.ITB',STATUS='OLD',FORM='UNFORMATTED')
               CALL SMSREAD(VINT,VINT2)
               CLOSE(51)
            ELSE
               IF (DOIT_TB .EQ. 1) THEN
                  OPEN(UNIT=51,FILE = NAME(1:J-1)//'.ITB',STATUS='UNKNOWN',FORM='UNFORMATTED')
               ENDIF
               CALL SMSNEW(DOIT_TB,VINT,VINT2)
            ENDIF
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
         GO TO 800
      ENDIF

      IF (AVAIL_OB) THEN
         J = INDEX(NAME,' ')
         OPEN(UNIT=50,FILE = NAME(1:J-1)//'.IOB',STATUS='OLD',FORM='UNFORMATTED')
         IF (YES2) THEN
            CALL DENSREAD_SELTZ_MPI(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6, &
                 DINT7,DINT1VEC,DENS1VEC,NRNUC)
         ELSE
            CALL DENSREAD_MPI(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
         END IF
         CLOSE(50)
      ELSE
         IF(YES2) THEN
            CALL DENSNEW_SELTZ_MPI(DOIT_OB,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6, &
                 DINT7,DINT1VEC,DENS1VEC,NRNUC)
         ELSE
            CALL DENSNEW_MPI(DOIT_OB,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
         END IF
      ENDIF

      IF (AVAIL_TB) THEN
         J = INDEX(NAME,' ')
         OPEN(UNIT=51,FILE = NAME(1:J-1)//'.ITB',STATUS='OLD',FORM='UNFORMATTED')
         CALL SMSREAD_MPI(VINT,VINT2)
         CLOSE(51)
      ELSE
         CALL SMSNEW_MPI(DOIT_TB,VINT,VINT2)
      ENDIF

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      CALL GDRSUMMPI_ROOT(SMSC1(1), NVEC)
      CALL GDRSUMMPI_ROOT(SMSC2(1), NVEC)
      CALL GDRSUMMPI_ROOT(DENS1(1), NVEC)
      CALL GDRSUMMPI_ROOT(DENS2(1), NVEC)
      CALL GDRSUMMPI_ROOT(DENS3(1), NVEC)
      CALL GDRSUMMPI_ROOT(DENS4(1), NVEC)
      CALL GDRSUMMPI_ROOT(DENS5(1), NVEC)
      CALL GDRSUMMPI_ROOT(DENS6(1), NVEC)
      CALL GDRSUMMPI_ROOT(DENS7(1), NVEC)
      IF (YES2) CALL GDRSUMMPI_ROOT(DENS1VEC(1,1), NVEC*NRNUC)

  800 CONTINUE
      IF (MYID .NE. 0) GO TO 900

      WRITE (24,'(a24,i3)') 'Number of eigenvalues: ',NVEC
      WRITE(24,*)
      WRITE(24,*)
      WRITE (24,309)
      DO I = 1, NVEC
         WRITE (24,323) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
              EAV+EVAL(I)
      END DO
      WRITE (24,308)
      DO I = 1,NVEC
         WRITE(24,314)
         WRITE (24,313) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
              DENS7(I),(DENS2(I)-DENS7(I)),DENS2(I)
         HzNMSu=DENS7(I)*AUMAMU*AUCM*CCMS
         WRITE (24,316) HzNMSu*(1.0D-09), &
              (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-09), &
              DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-09)
      END DO

      WRITE (24,302)
      DO I = 1,NVEC
         WRITE(24,314)
         WRITE (24,313) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
              SMSC1(I),SMSC2(I)-SMSC1(I),SMSC2(I)
         HzSMSu=SMSC1(I)*AUMAMU*AUCM*CCMS
         IF (DABS(HzSMSu) .LE. (10**8)) THEN
            WRITE (24,315) HzSMSu*(1.0D-06), &
                 (SMSC2(I)-SMSC1(I))*AUMAMU*AUCM*CCMS*(1.0D-06), &
                 SMSC2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
         ELSE
            WRITE (24,316) HzSMSu*(1.0D-09), &
                 (SMSC2(I)-SMSC1(I))*AUMAMU*AUCM*CCMS*(1.0D-09), &
                 SMSC2(I)*AUMAMU*AUCM*CCMS*(1.0D-09)
         ENDIF
      END DO
      IF(YES2) THEN
         WRITE (24,307)
         DO I = 1,NVEC
            CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT(:,I),RHO(I),RES(I))
            WRITE (24,303) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),RHO(I)
         END DO
         WRITE (24,317)
         DO I = 1,NVEC
            WRITE (24,324) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),FMAT(1,I),FMAT(2,I), &
                 FMAT(3,I),FMAT(4,I),RES(I)
         END DO
         WRITE (24,327)
         DO I = 1,NVEC
            WRITE (24,334) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),FMAT(5,I),FMAT(6,I)
         END DO
         WRITE(24,*)
      ELSE
         WRITE (24,306)
         DO I = 1,NVEC
            WRITE (24,303) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),DENS1(I)
         END DO
      END IF

  900 CONTINUE
      CALL DALLOC (SMSC1,'SMSC1','RIS_CAL_MPI')
      CALL DALLOC (SMSC2,'SMSC2','RIS_CAL_MPI')
      CALL DALLOC (DENS1,'DENS1','RIS_CAL_MPI')
      CALL DALLOC (DENS2,'DENS2','RIS_CAL_MPI')
      CALL DALLOC (DENS3,'DENS3','RIS_CAL_MPI')
      CALL DALLOC (DENS4,'DENS4','RIS_CAL_MPI')
      CALL DALLOC (DENS5,'DENS5','RIS_CAL_MPI')
      CALL DALLOC (DENS6,'DENS6','RIS_CAL_MPI')
      CALL DALLOC (DENS7,'DENS7','RIS_CAL_MPI')

      IF (MYID .EQ. 0) THEN
         WRITE(*,*) '-------------------------------'
         WRITE(*,*) 'RIS_CAL_MPI: Execution Finished ...'
         WRITE(*,*) '-------------------------------'
      ENDIF
      CALL STOPTIME (ncount1, 'RIS_CAL_MPI')
      RETURN
301   FORMAT (//' CUTOFF set to ',1PD22.15)
302   FORMAT (//' Level  J Parity  Specific mass shift parameter')
303   FORMAT (1X,I3,5X,2A4,3X,D20.10)
304   FORMAT (1X,I3,5X,2A4,3X,2D20.10)
324   FORMAT (1X,I3,5X,2A4,3X,4D20.10,F10.4)
334   FORMAT (1X,I3,5X,2A4,3X,2D20.10)
306   FORMAT (//' Electron density in atomic units' &
           //' Level  J Parity',8X,'DENS (a.u.)'/)
307   FORMAT (//' Level  J Parity  Electron density in atomic units' &
           //24X,'Dens. (a.u.)')
317   FORMAT (//' Level  J Parity  Field shift electronic factors and av&
           erage point discrepancy in fit' &
           //24X, &
           'F0 (GHz/fm^2)',6X,' F2 (GHz/fm^4)', &
           6X,' F4 (GHz/fm^6)',6X,' F6 (GHz/fm^8)' &
           6X,' Disc. (per mille)')
327   FORMAT (//' Level  J Parity  Field shift electronic factors (corre&
           cted for varying density inside nucleus)' &
           //24X, &
           'F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)')
308   FORMAT (//' Level  J Parity  Normal mass shift parameter')
309   FORMAT (' Level  J Parity  Energy')
313   FORMAT (1X,I3,5X,2A4,3x,3D20.10,'  (a.u.)')
333   FORMAT (1X,I3,5X,2A4,3x,3D20.10,'  (GHz u)')
323   FORMAT (1X,I3,5X,2A4,3x,D20.10,'  (a.u.)')
314   FORMAT (/,29X,'<K^1>',13X,'<K^2+K^3>',9X,'<K^1+K^2+K^3>')
315   FORMAT (20X,3D20.10 ,2X,'(MHz u)')
316   FORMAT (20X,3D20.10,2X,'(GHz u)')
      RETURN
      END SUBROUTINE RIS_CAL_MPI
