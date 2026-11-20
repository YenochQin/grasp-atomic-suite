      module setmcp_i
      implicit none (type, external)
      interface
      subroutine setmcp (core_subshell_count, input_symmetry_block_count, block_labels, filehead)
      implicit none (type, external)
      integer :: core_subshell_count
      integer, intent(in) :: input_symmetry_block_count
      character (len = 8) :: block_labels(*)
      character (len = *), intent(in) :: filehead
!...This routine performs I/O.
      end subroutine setmcp
      end interface
      end module setmcp_i
