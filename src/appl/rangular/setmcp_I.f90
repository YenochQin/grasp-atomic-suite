      module setmcp_i
      implicit none (type, external)
      interface
      subroutine setmcp (myid, nprocs, ncore, idblk, filehead)
      implicit none (type, external)
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      integer, intent(in) :: ncore
      character (len = 8), intent(in) :: idblk(*)
      character (len = *), intent(in) :: filehead
!...This routine performs I/O.
      end subroutine
      end interface
      end module
