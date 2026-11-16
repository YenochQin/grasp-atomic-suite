      module setcof_i
      implicit none (type, external)
      interface
      subroutine setcof (eol, j)
      implicit none (type, external)
      logical :: eol
      integer, intent(in) :: j
      end subroutine
      end interface
      end module
