      MODULE cteilsrk_C
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER, DIMENSION(:), pointer :: indteirk
      REAL(DOUBLE), DIMENSION(:), pointer :: valteirk

      CONTAINS

! GENINTRK numbers the cache entries from 1 and assigns entry i to
! MOD(i,nprocs). Count the local evaluations, not the replicated cache.
      PURE INTEGER FUNCTION RK_LOCAL_COUNT(N, MYID, NPROCS) RESULT(NLOCAL)
      INTEGER, INTENT(IN) :: N, MYID, NPROCS
      IF (N < 0 .OR. NPROCS < 1) ERROR STOP 'RK_LOCAL_COUNT: invalid rank layout'
      IF (MYID < 0 .OR. MYID >= NPROCS) ERROR STOP 'RK_LOCAL_COUNT: invalid rank layout'
      NLOCAL = N/NPROCS
      IF (MYID > 0 .AND. MYID <= MOD(N,NPROCS)) NLOCAL = NLOCAL + 1
      END FUNCTION RK_LOCAL_COUNT
      END MODULE cteilsrk_C
