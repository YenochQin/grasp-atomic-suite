      module estim_i
      implicit none (type, external)
      interface
      subroutine estim (j)
      implicit none (type, external)
      integer, intent(in) :: j
      end subroutine estim
      end interface
      end module estim_i
