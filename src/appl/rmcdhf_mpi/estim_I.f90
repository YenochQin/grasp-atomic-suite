      module estim_i
      implicit none (type, external)
      interface
      subroutine estim (j)
      implicit none (type, external)
      integer, intent(in) :: j
      end subroutine
      end interface
      end module
