      module genmat_i
      implicit none (type, external)
      interface
      subroutine genmat (atwinv, jblock, myid, nprocs, elsto, irestart, slf_en)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64) :: atwinv
!VAST...Dummy argument ATWINV is not referenced in this routine.
      integer, intent(in) :: jblock
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      real(real64), intent(out) :: elsto
      integer, intent(out) :: irestart
      real(real64) :: slf_en(*)
      end subroutine
      end interface
      end module
