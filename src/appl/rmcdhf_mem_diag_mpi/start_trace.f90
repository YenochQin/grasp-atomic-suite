! Read-only START diagnostics, private to rmcdhf_mem_diag_mpi.
      MODULE start_trace
      USE iso_fortran_env, ONLY: error_unit
      USE vast_kind_param, ONLY: DOUBLE
      USE p0_trace, ONLY: p0_trace_context
      USE orb_C, ONLY: NP, NH, NAK, E, GAMA
      USE grid_C, ONLY: R, RP, H
      USE def_C, ONLY: C, Z, ACCY
      USE npar_C, ONLY: NPARM
      USE pote_C, ONLY: XP, XQ, YP
      USE scf_C, ONLY: NDCOF, NDA, DA
      USE wave_C, ONLY: PZ, PF, QF
      USE cnc_C, ONLY: CNC6C
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: start_trace_open, start_trace_close, start_trace_snapshot
      LOGICAL :: enabled = .FALSE.
      INTEGER :: points_unit = 0, terms_unit = 0, rank_id = 0
      CONTAINS

      SUBROUTINE start_trace_open(directory, rank)
      CHARACTER(LEN=*), INTENT(IN) :: directory
      INTEGER, INTENT(IN) :: rank
      CHARACTER(LEN=16) :: option, rank_string
      INTEGER :: length, status, ios
      CALL start_trace_close()
      rank_id = rank
      option = '0'
      CALL GET_ENVIRONMENT_VARIABLE('GRASP_START_TRACE', LENGTH=length, STATUS=status)
      IF (status == 0 .AND. length > 0) THEN
         IF (length > LEN(option)) THEN
            CALL start_trace_error('invalid GRASP_START_TRACE')
            RETURN
         ENDIF
         CALL GET_ENVIRONMENT_VARIABLE('GRASP_START_TRACE', option)
      ENDIF
      SELECT CASE (TRIM(option))
      CASE ('0')
         RETURN
      CASE ('1')
         IF (rank /= 0) RETURN
      CASE ('all')
         CONTINUE
      CASE DEFAULT
         CALL start_trace_error('expected GRASP_START_TRACE=0, 1, or all')
         RETURN
      END SELECT
      WRITE(rank_string,'(I0)') rank
      OPEN(NEWUNIT=points_unit, FILE=TRIM(directory)//'/rmcdhf.start.rank'//TRIM(rank_string)//'.csv', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
         points_unit = 0
         CALL start_trace_error('cannot open START points CSV')
         RETURN
      ENDIF
      OPEN(NEWUNIT=terms_unit, FILE=TRIM(directory)//'/rmcdhf.start_terms.rank'//TRIM(rank_string)//'.csv', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
         terms_unit = 0
         CALL start_trace_error('cannot open START coefficients CSV')
         RETURN
      ENDIF
      enabled = .TRUE.
      WRITE(points_unit,'(A)',IOSTAT=ios) &
         'cycle,update,rank,orbital,label,solve_iteration,itype,stage,start_iteration,point,'// &
         'r,rp,h,c,z,nparm,p0,q0,gamma,kappa,energy,accy,ndcof,xp,xq,yp,rsep,rseq,rdp,rdq,'// &
         'p_scaled,q_scaled,stored_p,stored_q,delta_scaled,delta_physical,meets_accy,fallback'
      IF (.NOT.start_trace_io_ok(ios)) RETURN
      WRITE(terms_unit,'(A)',IOSTAT=ios) &
         'cycle,update,rank,orbital,label,solve_iteration,itype,kind,index,peer,value,peer_p0'
      IF (.NOT.start_trace_io_ok(ios)) RETURN
      ! This literal also lets the job preflight reject an older executable.
      IF (rank == 0) WRITE(*,'(A)') 'GRASP_START_TRACE_V1 enabled'
      END SUBROUTINE start_trace_open

      SUBROUTINE start_trace_close()
      INTEGER :: ios
      IF (points_unit /= 0) CLOSE(points_unit, IOSTAT=ios)
      IF (terms_unit /= 0) CLOSE(terms_unit, IOSTAT=ios)
      points_unit = 0
      terms_unit = 0
      enabled = .FALSE.
      END SUBROUTINE start_trace_close

      SUBROUTINE start_trace_error(message)
      CHARACTER(LEN=*), INTENT(IN) :: message
      WRITE(error_unit,'(A)') '[START_TRACE] '//message//'; START tracing disabled'
      CALL start_trace_close()
      END SUBROUTINE start_trace_error

      LOGICAL FUNCTION start_trace_io_ok(status) RESULT(ok)
      INTEGER, INTENT(IN) :: status
      INTEGER :: ios
      ios = status
      IF (ios == 0) FLUSH(points_unit, IOSTAT=ios)
      IF (ios == 0) FLUSH(terms_unit, IOSTAT=ios)
      ok = ios == 0
      IF (.NOT.ok) CALL start_trace_error('output error')
      END FUNCTION start_trace_io_ok

      SUBROUTINE start_trace_snapshot(stage, orbital, itype, iteration, pzero, qzero, p, q, &
                                     rsep, rseq, rdp, rdq, raw_delta, delta, fallback)
      CHARACTER(LEN=*), INTENT(IN) :: stage
      INTEGER, INTENT(IN) :: orbital, itype, iteration
      REAL(DOUBLE), INTENT(IN) :: pzero, qzero, p(:), q(:), rsep(6), rseq(6), rdp(6), rdq(6)
      REAL(DOUBLE), INTENT(IN) :: raw_delta, delta
      LOGICAL, INTENT(IN) :: fallback
      INTEGER :: cycle, update, solve_iteration, i, j, peer, ios
      CHARACTER(LEN=16) :: label
      IF (.NOT.enabled) RETURN
      CALL p0_trace_context(cycle, update, solve_iteration)
      WRITE(label,'(I0,A)') NP(orbital), TRIM(NH(orbital))
      DO i = 1, 6
         WRITE(points_unit,'(*(G0,:,","))',IOSTAT=ios) &
            cycle, update, rank_id, orbital, TRIM(label), solve_iteration, itype, stage, iteration, i, &
            R(i), RP(i), H, C, Z, NPARM, pzero, qzero, GAMA(orbital), NAK(orbital), E(orbital), ACCY, &
            NDCOF, XP(i), XQ(i), YP(i), rsep(i), rseq(i), rdp(i), rdq(i), p(i), q(i), &
            PF(i,orbital), QF(i,orbital), raw_delta, delta, iteration > 0 .AND. delta <= ACCY, fallback
         IF (.NOT.start_trace_io_ok(ios)) RETURN
      END DO
      IF (stage /= 'initial') RETURN
      IF (itype == 2) THEN
         DO i = 1, NDCOF
            peer = NDA(i)
            WRITE(terms_unit,'(*(G0,:,","))',IOSTAT=ios) &
               cycle, update, rank_id, orbital, TRIM(label), solve_iteration, itype, &
               'coupling', i, peer, DA(i), PZ(peer)
            IF (.NOT.start_trace_io_ok(ios)) RETURN
         END DO
      ENDIF
      DO j = 2, 6
         DO i = 1, 6
            WRITE(terms_unit,'(*(G0,:,","))',IOSTAT=ios) &
               cycle, update, rank_id, orbital, TRIM(label), solve_iteration, itype, &
               'cnc', i, j, CNC6C(i,j), 0.0D0
            IF (.NOT.start_trace_io_ok(ios)) RETURN
         END DO
      END DO
      END SUBROUTINE start_trace_snapshot
      END MODULE start_trace
