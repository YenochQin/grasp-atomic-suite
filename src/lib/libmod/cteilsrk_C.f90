      module cteilsrk_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(:), pointer :: indteirk
      real(double), dimension(:), pointer :: valteirk

      contains

! GENINTRK numbers the cache entries from 1 and assigns entry i to
! MOD(i,nprocs). Count the local evaluations, not the replicated cache.
      pure integer function rk_local_count(n, myid, nprocs) result(nlocal)
      integer, intent(in) :: n, myid, nprocs
      if (n < 0 .or. nprocs < 1) error stop 'RK_LOCAL_COUNT: invalid rank layout'
      if (myid < 0 .or. myid >= nprocs) error stop 'RK_LOCAL_COUNT: invalid rank layout'
      nlocal = n/nprocs
      if (myid > 0 .and. myid <= mod(n,nprocs)) nlocal = nlocal + 1
      end function rk_local_count
      end module cteilsrk_c
