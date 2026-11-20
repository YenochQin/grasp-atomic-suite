      module setcsl_i
      implicit none (type, external)
      interface
      subroutine setcsl (name, core_subshell_count, input_symmetry_block_count, block_labels)
      implicit none (type, external)
      character (len = *) :: name
      integer :: core_subshell_count
      integer, intent(in) :: input_symmetry_block_count
      character (len = 8) :: block_labels(*)
      end subroutine setcsl
      end interface
      end module setcsl_i
