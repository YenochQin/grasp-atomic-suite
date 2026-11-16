      module setham_i
      implicit none (type, external)
      interface
      subroutine setham (myid, nprocs, jblock, elsto,icstrt, nelmntt,  &
                         atwinv,slf_en)
      use, intrinsic :: iso_fortran_env, only: real64, int64
      implicit none (type, external)
      integer(int64)       :: nelmntt
      integer             :: jblock, icstrt
      integer, intent(in) :: myid, nprocs
      real(real64) :: elsto, atwinv
      real(real64) :: slf_en(*)
      end subroutine setham
      end interface
      end module
