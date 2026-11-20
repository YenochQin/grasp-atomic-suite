      module consis_i
      implicit none (type, external)
      interface
      subroutine consis (orbital_index)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      end subroutine consis
      end interface
      end module consis_i
