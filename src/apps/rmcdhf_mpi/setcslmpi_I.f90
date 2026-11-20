      module setcslmpi_i
      implicit none (type, external)
      interface
      subroutine setcslmpi (name, core_subshell_count, block_labels)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      integer :: core_subshell_count
      character (len = 8) :: block_labels(*)
!VAST...Dummy argument IDBLK is not referenced in this routine.
      end subroutine setcslmpi
      end interface
      end module setcslmpi_i
