      module hmout_i
      implicit none (type, external)
      interface
      subroutine hmout (restart_matrix_unit)
      implicit none (type, external)
      integer, intent(in) :: restart_matrix_unit
      end subroutine hmout
      end interface
      end module hmout_i
