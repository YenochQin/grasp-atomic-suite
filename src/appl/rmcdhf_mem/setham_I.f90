      module setham_i
      implicit none (type, external)
      interface
      subroutine setham (jblock, myid, nprocs, ncontr_tot, ncontr_tot2)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      integer, intent(in) :: jblock
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      integer(int64) :: ncontr_tot(20), ncontr_tot2(20)
      end subroutine setham
      end interface
      end module setham_i
