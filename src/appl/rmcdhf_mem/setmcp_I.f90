      module setmcp_i
      implicit none (type, external)
      interface
      subroutine setmcp (ncore, nblkin, idblk, filehead)
      implicit none (type, external)
      integer :: ncore
      integer, intent(in) :: nblkin
      character (len = 8) :: idblk(*)
      character (len = *), intent(in) :: filehead
!...This routine performs I/O.
      end subroutine setmcp
      end interface
      end module setmcp_i
