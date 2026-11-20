      module matrix_i
      implicit none (type, external)
      interface
      subroutine matrix (unused_core_subshell_count, max_orbital_two_j)
      implicit none (type, external)
      integer,  intent(in):: unused_core_subshell_count, max_orbital_two_j
      end subroutine matrix
      end interface
      end module matrix_i
