      module setmcp2_i
      implicit none (type, external)
      interface
      subroutine setmcp2 (myid, nprocs, ncore, idblk, filehead)
      implicit none (type, external)
      integer, intent(in) :: myid
      integer :: nprocs
      integer :: ncore
      character (len = 8) :: idblk(*)
      character (len = *) :: filehead
      end subroutine
      end interface
      end module
