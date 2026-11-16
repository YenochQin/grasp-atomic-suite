      module defcor_i
      implicit none (type, external)
      interface
      subroutine defcor (j)
      implicit none (type, external)
      integer, intent(in) :: j
      end subroutine
      end interface
      end module
