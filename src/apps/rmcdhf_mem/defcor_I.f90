      module defcor_i
      implicit none (type, external)
      interface
      subroutine defcor (orbital_index)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      end subroutine defcor
      end interface
      end module defcor_i
