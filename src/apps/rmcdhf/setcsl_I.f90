      module setcsl_i
      implicit none (type, external)
      interface
      subroutine setcsl (name, core_subshell_count, block_labels)
      implicit none (type, external)
      character (len = *), intent(in) :: name
      integer :: core_subshell_count
      character (len = 8) :: block_labels(*)
      end subroutine setcsl
      end interface
      end module setcsl_i
