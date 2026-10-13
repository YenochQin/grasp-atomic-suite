! Read-only tracing; no changes to numerical state or MPI collectives.
      MODULE p0_trace
      USE iso_fortran_env, ONLY: error_unit
      USE vast_kind_param, ONLY: DOUBLE
      USE orb_C, ONLY: NP, NH, NAK, E, GAMA
      USE wave_C, ONLY: PZ, PF, QF, MF
      USE grid_C, ONLY: R
      USE scf_C, ONLY: METHOD, SCNSTY
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: p0_trace_open, p0_trace_close, p0_trace_cycle
      PUBLIC :: p0_trace_update, p0_trace_iteration
      PUBLIC :: p0_trace_stored, p0_trace_trial
      LOGICAL :: enabled = .FALSE.
      INTEGER :: trace_unit = 0, rank_id = 0, cycle_id = 0
      INTEGER :: update_id = 0, active_orbital = 0, inner_iteration = 0
      CONTAINS

      SUBROUTINE p0_trace_open(directory, rank)
      CHARACTER(LEN=*), INTENT(IN) :: directory
      INTEGER, INTENT(IN) :: rank
      CHARACTER(LEN=16) :: option, rank_string
      CHARACTER(LEN=:), ALLOCATABLE :: filename
      INTEGER :: status, length, ios
      CALL p0_trace_close()
      rank_id = rank
      cycle_id = 0
      update_id = 0
      active_orbital = 0
      inner_iteration = 0
      option = '1'
      CALL GET_ENVIRONMENT_VARIABLE('GRASP_P0_TRACE', LENGTH=length, STATUS=status)
      IF (status == 0 .AND. length > 0) THEN
         IF (length > LEN(option)) THEN
            IF (rank == 0) WRITE(error_unit,'(A)') '[P0_TRACE] invalid GRASP_P0_TRACE; tracing disabled'
            RETURN
         ENDIF
         CALL GET_ENVIRONMENT_VARIABLE('GRASP_P0_TRACE', option)
      ENDIF
      SELECT CASE (TRIM(option))
      CASE ('0')
         RETURN
      CASE ('1')
         IF (rank /= 0) RETURN
      CASE ('all')
         CONTINUE
      CASE DEFAULT
         IF (rank == 0) WRITE(error_unit,'(A)') '[P0_TRACE] expected GRASP_P0_TRACE=0, 1, or all'
         RETURN
      END SELECT
      WRITE(rank_string,'(I0)') rank
      filename = TRIM(directory)//'/rmcdhf.p0.rank'//TRIM(rank_string)//'.csv'
      OPEN(NEWUNIT=trace_unit, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(error_unit,'(A)') '[P0_TRACE] cannot open '//filename//'; tracing disabled'
         RETURN
      ENDIF
      enabled = .TRUE.
      WRITE(trace_unit,'(A)',IOSTAT=ios) &
         'cycle,update,rank,stage,active_orbital,orbital,label,n,kappa,inner_iteration,peer,method,mtp,'// &
         'p0,energy,gamma,scnsty,overlap,scale,norm,mixing,alpha,match_denominator,'// &
         'p2,p3,p4,p5,p6,p7,q2,q3,q4,q5,q6,q7,r2,r3,r4,r5,r6,r7'
      CALL trace_io_check(ios)
      END SUBROUTINE p0_trace_open

      SUBROUTINE p0_trace_close()
      INTEGER :: ios
      IF (enabled) CLOSE(trace_unit, IOSTAT=ios)
      enabled = .FALSE.
      trace_unit = 0
      END SUBROUTINE p0_trace_close

      SUBROUTINE p0_trace_cycle(cycle)
      INTEGER, INTENT(IN) :: cycle
      cycle_id = cycle
      update_id = 0
      active_orbital = 0
      inner_iteration = 0
      END SUBROUTINE p0_trace_cycle

      SUBROUTINE p0_trace_update(orbital)
      INTEGER, INTENT(IN) :: orbital
      update_id = update_id + 1
      active_orbital = orbital
      inner_iteration = 0
      CALL p0_trace_stored('update_begin', orbital)
      END SUBROUTINE p0_trace_update

      SUBROUTINE p0_trace_iteration(iteration)
      INTEGER, INTENT(IN) :: iteration
      inner_iteration = iteration
      END SUBROUTINE p0_trace_iteration

      SUBROUTINE p0_trace_stored(stage, orbital, peer, overlap, scale, norm, mixing)
      CHARACTER(LEN=*), INTENT(IN) :: stage
      INTEGER, INTENT(IN) :: orbital
      INTEGER, OPTIONAL, INTENT(IN) :: peer
      REAL(DOUBLE), OPTIONAL, INTENT(IN) :: overlap, scale, norm, mixing
      IF (.NOT.enabled) RETURN
      CALL p0_trace_trial(stage, orbital, PZ(orbital), PF(:,orbital), QF(:,orbital), &
         MF(orbital), peer=peer, overlap=overlap, scale=scale, norm=norm, mixing=mixing)
      END SUBROUTINE p0_trace_stored

      SUBROUTINE p0_trace_trial(stage, orbital, pzero, p, q, points, peer, &
                               overlap, scale, norm, mixing, alpha, match_denominator)
      CHARACTER(LEN=*), INTENT(IN) :: stage
      INTEGER, INTENT(IN) :: orbital, points
      REAL(DOUBLE), INTENT(IN) :: pzero, p(:), q(:)
      INTEGER, OPTIONAL, INTENT(IN) :: peer
      REAL(DOUBLE), OPTIONAL, INTENT(IN) :: overlap, scale, norm, mixing, alpha, match_denominator
      REAL(DOUBLE) :: ps(6), qs(6), details(6)
      INTEGER :: k, other, ios
      CHARACTER(LEN=16) :: label
      IF (.NOT.enabled) RETURN
      ps = 0.0D0
      qs = 0.0D0
      DO k = 2, MIN(7, points, SIZE(p))
         ps(k-1) = p(k)
      END DO
      DO k = 2, MIN(7, points, SIZE(q))
         qs(k-1) = q(k)
      END DO
      details = 0.0D0
      IF (PRESENT(overlap)) details(1) = overlap
      IF (PRESENT(scale)) details(2) = scale
      IF (PRESENT(norm)) details(3) = norm
      IF (PRESENT(mixing)) details(4) = mixing
      IF (PRESENT(alpha)) details(5) = alpha
      IF (PRESENT(match_denominator)) details(6) = match_denominator
      other = 0
      IF (PRESENT(peer)) other = peer
      WRITE(label,'(I0,A)') NP(orbital), TRIM(NH(orbital))
      WRITE(trace_unit,'(*(G0,:,","))',IOSTAT=ios) cycle_id, update_id, rank_id, TRIM(stage), &
         active_orbital, orbital, TRIM(label), NP(orbital), NAK(orbital), inner_iteration, other, &
         METHOD(orbital), points, pzero, E(orbital), GAMA(orbital), SCNSTY(orbital), &
         details, ps, qs, R(2:7)
      CALL trace_io_check(ios)
      END SUBROUTINE p0_trace_trial

      SUBROUTINE trace_io_check(status)
      INTEGER, INTENT(IN) :: status
      INTEGER :: ios
      IF (.NOT.enabled) RETURN
      ios = status
      IF (ios == 0) FLUSH(trace_unit, IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE(error_unit,'(A)') '[P0_TRACE] output error; tracing disabled for this rank'
         CALL p0_trace_close()
      ENDIF
      END SUBROUTINE trace_io_check
      END MODULE p0_trace
