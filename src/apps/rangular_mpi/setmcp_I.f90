      module setmcp_i
      implicit none (type, external)
      interface
!cjb  myid, nprocs = NOT args
      subroutine setmcp (core_subshell_count, block_labels, filehead)
      implicit none (type, external)
!     INTEGER, INTENT(IN) :: MYID
!     INTEGER, INTENT(IN) :: NPROCS
      integer, intent(in) :: core_subshell_count
      character (len = 8), intent(in) :: block_labels(*)
      character (len = *), intent(in) :: filehead
!...This routine performs I/O.
      end subroutine setmcp
      end interface
      end module setmcp_i
