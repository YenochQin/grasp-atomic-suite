      module estim_i
      implicit none (type, external)
      interface
      subroutine estim (orbital_index)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      end subroutine estim
      end interface
      end module estim_i
